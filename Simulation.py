import numpy as np
import matplotlib.pyplot as plt

class RInSimulation:
    def __init__(self, h, stim_amp=-1.0, stim_delay=100.0, stim_dur=800.0, tstop=1000.0):
        """
        Initialize the simulation with cell and stimulation parameters.
        
        Parameters:
          cell     : the NEURON cell object (expected to have a 'soma' section)
          stim_amp : amplitude of the injected current (nA)
          stim_delay: delay before the stimulus begins (ms)
          stim_dur  : duration of the stimulus (ms)
          tstop     : simulation end time (ms)
          dt        : simulation time step (ms)
        """
        self.h = h
        self.h.tstop = tstop
        self.stim_amp = stim_amp
        self.stim_delay = stim_delay
        self.stim_dur = stim_dur
        self.stim = None      # will hold the IClamp object
        self.t_vec = None     # time recorder
        self.v_vec = None     # voltage recorder

    def setup_stimulation(self):
        """Set up the IClamp at the midpoint of the soma."""
        self.stim = self.h.IClamp(self.h.soma[0](0.5))
        self.stim.amp = self.stim_amp
        self.stim.delay = self.stim_delay
        self.stim.dur = self.stim_dur

    def setup_recording(self):
        """Set up vectors to record time and somatic voltage."""
        self.t_vec = self.h.Vector()
        self.v_vec = self.h.Vector()
        self.t_vec.record(self.h._ref_t)
        self.v_vec.record(self.h.soma[0](0.5)._ref_v)

    def run_simulation(self):
        """Initialize and run the simulation."""
        self.setup_stimulation()
        self.setup_recording()
        # You may need to adjust the initial voltage (here, set to -65 mV)
        self.h.finitialize(-65)
        self.h.run()

    def plot_voltage(self):
        """Plot the recorded voltage trace."""
        t = np.array(self.t_vec)
        v = np.array(self.v_vec)
        plt.figure()
        plt.plot(t, v)
        plt.xlabel("Time (ms)")
        plt.ylabel("Voltage (mV)")
        plt.title("Voltage Trace")
        plt.savefig("voltage_trace_Rin.png")

    def measure_r_in(self):
        """
        Run the simulation, plot the voltage trace, and calculate the input resistance.
        
        The calculation uses:
          R_in = (V_rest - V_trough) / abs(I_stim)
        where:
          - V_rest is the average voltage before stimulation (over a 10 ms window).
          - V_trough is the minimum voltage during the stimulus.
        """
        self.run_simulation()
        self.plot_voltage()
        t = np.array(self.t_vec)
        v = np.array(self.v_vec)
        dt = self.h.dt

        # Determine the indices corresponding to the start and end of the stimulus
        stim_start_idx = int(self.stim_delay / dt)
        stim_end_idx = int((self.stim_delay + self.stim_dur) / dt)

        # determine indices to use for V_rest and V_trough
        V_rest_idx = stim_start_idx - 1
        V_trough_idx = stim_end_idx - 1

        # find the voltages
        V_rest = v[V_rest_idx]
        V_trough = v[V_trough_idx]
      
        # find the actual time
        V_rest_time = V_rest_idx * dt
        V_trough_time = V_trough_idx * dt
        
        # Calculate input resistance in MOhm (mV/nA)
        r_in = (V_rest - V_trough) / self.stim_amp

        print(f"V_rest [ {V_rest:.3} ] mV occurs at [ {V_rest_time:.3} ] ms")
        print(f"V_trough [ {V_trough:.3} ] mV occurs at [ {V_trough_time:.3} ] ms")
        print(f"r_in = [ {V_rest:.3} - {V_trough:.3} ] mV / [ {self.stim_amp:.3} ] nA  = {r_in:.3} MOhm")
      
        return r_in
