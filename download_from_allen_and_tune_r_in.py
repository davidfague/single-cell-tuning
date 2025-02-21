# example use: python download_from_allen_and_tune_r_in 488683425 "user_specifications.json"
# or: python download_from_allen_and_tune_r_in "http://celltypes.brain-map.org/experiment/electrophysiology/488683425" "user_specifications.json"
# or python download_from_allen_and_tune_r_in 488683425
# or python download_from_allen_and_tune_r_in "http://celltypes.brain-map.org/experiment/electrophysiology/488683425"



from allensdk.api.queries.biophysical_api import BiophysicalApi
from allensdk.model.biophys_sim.config import Config
from allensdk.model.biophysical.utils import Utils
import re
import json
import numpy as np
import matplotlib.pyplot as plt
import subprocess

from Simulation import RInSimulation

def robust_int_conversion(input_value):
  """
  Robustly converts various input types to an integer.

  Args:
    input_value: The input value to convert.  Can be a string (e.g., URL,
      "488683425"), an integer, or other types.

  Returns:
     An integer representation of the input.

  Raises:
    ValueError: If the input cannot be converted to an integer.
    TypeError: If the input type is not supported.
  """
  if isinstance(input_value, int):
    return input_value
  elif isinstance(input_value, str):
    # Try to find a sequence of digits using a regular expression
    match = re.search(r'\d+', input_value)
    if match:
      try:
        return int(match.group(0))
      except ValueError:
        raise ValueError(f"Invalid digit sequence found in string: {input_value}")
    else:
      raise ValueError(f"No digits found in string: {input_value}")
  else:
    try:
      return int(input_value)
    except (ValueError, TypeError) as e:
      raise TypeError(f"Unsupported input type or invalid value: {input_value}") from e

def load_dictionary_from_json(file_path):
  """Loads a dictionary from a JSON file.

  Args:
    file_path: The path to the JSON file.

  Returns:
    A dictionary loaded from the JSON file.

  Raises:
    ValueError: If the JSON file is invalid or cannot be parsed.
    FileNotFoundError: If the specified file does not exist.
  """
  try:
    with open(file_path, 'r') as file:
      try:
        data = json.load(file)
        if not isinstance(data, dict):
          raise ValueError("JSON file does not contain a dictionary.")
        return data
      except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON format: {e}") from e
  except FileNotFoundError:
    raise FileNotFoundError(f"File not found: {file_path}")

def update_missing_passive_values(utils, user_specifications):
  # update missing properties to user_specifications if they're not already in the allen specifications
  if "e_pas" not in utils.description.data["passive"][0].keys():
    utils.description.data["passive"][0]["e_pas"] = user_specifications["e_pas"]
  
  if "cm" not in utils.description.data["passive"][0].keys():
    utils.description.data["passive"][0]["cm"] = user_specifications["cm"]
  
  if "ra" not in utils.description.data["passive"][0].keys():
    utils.description.data["passive"][0]["ra"] = user_specifications["ra"]

  return utils

def measure_soma_surface_area(h):
  # calculate soma surface area
  # in NEURON each segment is cylindrical, however
  # in NEURON, does not include the flat surfaces in surface area calculation.
  
  # circumference times length
  soma_surface_area = h.soma[0].diam * np.pi * h.soma[0].L # L in microns, diam in microns
  soma_surface_area = soma_surface_area * 1e-8 # convert from microns^2 to cm^2
  print(f"soma_surface_area {soma_surface_area:.3} cm^2")
  return soma_surface_area

def estimate_gbar_leak_for_user_spec_rin(soma_surface_area, user_specs_dict):
  # This estimate is worse for multicompartment models.
  # This estimate is worse for models with active channels.
  # This estimate is only used to illustrate an example tuning workflow.
  
  gbar_leak_estimate = (1 / user_specs_dict['R-in']) / soma_surface_area
  # 1 / (MOhm * cm^2) = uS / cm^2
  print(f"gbar_leak_estimate {gbar_leak_estimate:.5} uS / cm2") # needs to be in uS / cm2
  return gbar_leak_estimate

def update_sections(value_to_assign: float, data: dict, sections: list = ['soma'], var_to_update: list = ['g_pas']):
  original_entries = []
  corresponding_new_entries = []

  original_assignments = data.copy()
  new_assignments = original_assignments.copy()

  print(original_assignments)

  for entry in new_assignments: # each entry is the assignment of a conductance somewhere
    if (entry['section'] in sections) or ('all' in sections):
      # print(f"entry: {entry}")
      if (entry['name'] in var_to_update) or ('all' in var_to_update):
        print(f" updating entry: {entry}")
        original_entries.append(entry.copy())

        entry['value'] = value_to_assign
        corresponding_new_entries.append(entry)

  for i,orig_entry in enumerate(original_entries):
    print(f"updating {entry['section']} {entry['name']} from {entry['value']:.3} to {corresponding_new_entries[i]['value']:.3}  percent change: {(((corresponding_new_entries[i]['value'] - entry['value']) / entry['value']) * 100):.3}")

  # return original_assignments, new_assignments, original_entries, corresponding_new_entries
  return new_assignments

if __name__ == "__main__":
  if len(sys.argv) > 1:
    int = sys.argv[1]
    print(f"getting cell {int}")
  else:
    raise(ValueError("No argument provided."))

  # load user specifications
  if len(sys.argv) > 2: # third argument expecting a json
    user_specs_dict = load_dictionary_from_json(sys.argv[2])
  else:
    user_specs_dict = None
    
  bp = BiophysicalApi()
  query = bp.get_neuronal_models(int)

  bp.cache_stimulus = False # Change to False to not download the large stimulus NWB file
  bp.cache_data(query[0]['id']) # 'id'

  # compile the downloaded modfiles
  subprocess.run("nrnivmodl modfiles", shell=True, check=True)

  # Create the h object
  description = Config().load('manifest.json')
  utils = Utils(description)
  h = utils.h

  # convert 'values' from string to float
  for dict in utils.description.data['genome']:
    for key,item in dict.items():
      if key == 'value':
        dict[key] = float(item)

  utils = update_missing_passive_values(utils, user_specifications)

  # read morphology
  manifest = description.manifest
  morphology_path = description.manifest.get_path('MORPHOLOGY')
  utils.generate_morphology(morphology_path.encode('ascii', 'ignore'))

  # build the cell. Its parts will be assigned to the h object
  utils.load_cell_parameters()

  soma_surface_ara = measure_soma_surface_area(h)

  gbar_leak_estimate = estimate_gbar_leak_for_user_spec_rin(soma_surface_area, user_specs_dict)

  utils.description.data["genome"] = update_sections(
                  value_to_assign = gbar_leak_estimate,
                  sections = ['soma'],
                  var_to_update = ['g_pas'], 
                  data = utils.description.data["genome"])
  
  utils.load_cell_parameters()

  r_in_sim_obj = RInSimulation(h)
  r_in = r_in_sim_obj.measure_r_in()
