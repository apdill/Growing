import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import io
import math
from tqdm import tqdm
from scipy.stats import powerlaw

class Neuron:
    def __init__(self, potential, position, size, color='black', reset_potential=-1):
        self.position = position
        self.size = size
        self.color = color
        self.reset_potential = reset_potential
        self.potential = potential
        self.spiking = np.random.rand() > 0.5
        self.firing_rate = int(self.spiking) 
        self.nspikes = int(self.spiking)
        if self.spiking: self.potential = reset_potential
        
    def determine_neuron_color(self):
            if self.spiking: 
                color = 'green' 
                self.nspikes += 1
            else: 
                voltage_normalized = (self.potential + 1) / 2
                color_rgb = plt.cm.winter(voltage_normalized)[:3]
                color = color_rgb
            return color
        
class Network:
    def __init__(self, N=60, tfr = 0.01, sandbox_size=(1000, 1000), size_range=(5, 50)):
        
        self.N = N
        self.sandbox_size = sandbox_size
        self.size_range = size_range
        self.neurons = self.init_neurons()
        
        self.M = self.init_connectivity_matrix()
        triu = np.triu_indices_from(self.M, k=1)
        tril = np.tril_indices_from(self.M, k=-1)
        
        self.indices = (np.concatenate((triu[0],tril[0])), np.concatenate((triu[1],tril[1])))
        
        self.update_connectivity_matrix()# evolving connectivity matrix M
        
        self.tfr = tfr # target firing rate
        self.mfr = 0.5
        
        self.growth_coeff = 0.01
        self.t = 0
      

    def init_connectivity_matrix(self):
        M = np.random.rand(self.N, self.N)
        M = M - np.transpose(M)
        np.fill_diagonal(M, 0)
        M = M / np.sum(abs(M))
        return M


    def init_neurons(self):
        sandbox_width, sandbox_height = self.sandbox_size
        neurons = []
        for _ in range(self.N):
            potential = 2*np.random.rand()-1
            position = tuple(np.random.rand(2) * [sandbox_width, sandbox_height])
            size = np.random.uniform(self.size_range[0], self.size_range[1])
            color = 'black'
            neurons.append(Neuron(potential, position, size, color))
        return neurons
          
            
    def check_overlap(self, neuron1, neuron2):
        distance_squared = (neuron1.position[0] - neuron2.position[0])**2 + (neuron1.position[1] - neuron2.position[1])**2
        return (neuron1.size + neuron2.size)**2 > distance_squared


    def circle_overlap_area(self, neuron1, neuron2):
        x1, y1, r1 = neuron1.position[0], neuron1.position[1], neuron1.size
        x2, y2, r2 = neuron2.position[0], neuron2.position[1], neuron2.size
        d = ((x2 - x1)**2 + (y2 - y1)**2)
        
        if d <= (r1 - r2)**2 and r1 >= r2:
            return math.pi * r2**2
        
        elif d <= (r1 - r2)**2 and r1 < r2:
            return math.pi * r1**2
        
        else:
            d = math.sqrt(d)
            h = (1/(2*d))*math.sqrt(2*(d**2)*(r1**2)+2*(d**2)*(r2**2) + 2*(r1**2)*(r2**2)-d**4 - r1**4 - r2**4)
            
            return (r1**2)*math.asin(h/r1)+(r2**2)*math.asin(h/r2) - d*h
        
    def update_connectivity_matrix(self):
        
        for k in range(len(self.indices[0])):
            i = self.indices[0][k]
            j = self.indices[1][k]
            
            if self.check_overlap(self.neurons[i], self.neurons[j]):
                overlap_area = self.circle_overlap_area(self.neurons[i], self.neurons[j])
                self.M[i][j] = overlap_area/(math.pi*(self.size_range[1])**2)
                #self.M[j][i] = - self.M[i][j]
            else: 
                self.M[i][j] = 0
                #self.M[j][i] = 0
                      
            


    def simulate(self, duration = 100, dt = 0.1, spike_input=1, input_noise = 0.1, leak=0.01, spike_threshold = 1, REST=-1, mfr = False, animate=False):
        
        # initializing global lists
        
        spikes_list = []
        V_list = []
        firing_rates_list = []
        size_list = []
        
        # initializing V and spike state vectors
        
        V = np.array([neuron.potential for neuron in self.neurons])
        spikes = [neuron.spiking for neuron in self.neurons]
    
        
        iterations = int(duration/dt) 
        
        if animate:
            fig, ax = plt.subplots()
            frames = []
        
            for _ in tqdm(range(iterations), desc="Simulating"):
                
                self.t += dt
                
                # Updating V (potential). Spike inputs are translated immediately
                
                spike_input_contribution = spike_input * np.dot(self.M, spikes)
                dV_dt = input_noise * np.random.rand(self.N) - leak * (V - REST)
                V += spike_input_contribution + (dV_dt * dt)
                
                # Checking which neurons spike, then resetting potentials
                
                spikes = (V >= spike_threshold)
                V[spikes] = REST
                V[V < REST] = REST

                # Updating neuron object attributes for neurons in network
                
                for i, neuron in enumerate(self.neurons):
                    neuron.spiking = spikes[i]
                    neuron.potential = V[i]
                    neuron.color = neuron.determine_neuron_color()
                    neuron.spiking = False
                    neuron.firing_rate = neuron.nspikes / self.t
                                    
                    if neuron.firing_rate < self.tfr:
                        neuron.size *= (1 + self.growth_coeff)
                    elif neuron.firing_rate >= self.tfr: 
                        neuron.size *= (1 - self.growth_coeff)
                    if neuron.size > self.size_range[1]: neuron.size = self.size_range[1]
                    if neuron.size < self.size_range[0]: neuron.size = self.size_range[0]
                
                # Updating global lists (keep track of all attributes for every iteration)
                
                self.update_connectivity_matrix()
                firing_rates = [neuron.firing_rate for neuron in self.neurons]
                neuron_sizes = [neuron.size for neuron in self.neurons]
                mean_firing_rate = np.mean(np.array(firing_rates))
                
                V_list.append(V.copy())
                spikes_list.append(spikes.copy())
                firing_rates_list.append(firing_rates.copy())
                size_list.append(neuron_sizes.copy())                
                
                ax.clear()
                for neuron in self.neurons:
                    circle = plt.Circle(neuron.position, neuron.size, color=neuron.color, alpha=0.5)
                    ax.add_artist(circle)
                ax.set_xlim(0, self.sandbox_size[0])
                ax.set_ylim(0, self.sandbox_size[1])
                ax.set_aspect('equal')
                ax.axis('off')
                buf = io.BytesIO()
                fig.savefig(buf, format='png')
                buf.seek(0)
                frames.append(imageio.imread(buf))

            imageio.mimsave('area_plasticity.gif', frames, duration=0.1)
            print('mean firing rate:', mean_firing_rate)    
            return spikes_list, V_list, firing_rates_list, size_list
        
        else:
            
            for _ in tqdm(range(iterations), desc="Simulating"):
                
                self.t += dt
                
                # Updating V (potential). Spike inputs are translated immediately
                
                spike_input_contribution = spike_input * np.dot(self.M, spikes)
                dV_dt = input_noise * np.random.rand(self.N) - leak * (V - REST)
                V += spike_input_contribution + (dV_dt * dt)
                
                # Checking which neurons spike, then resetting potentials
                
                spikes = (V >= spike_threshold)
                V[spikes] = REST
                V[V < REST] = REST

                # Updating neuron object attributes for neurons in network
                
                for i, neuron in enumerate(self.neurons):
                    neuron.spiking = spikes[i]
                    if neuron.spiking: neuron.nspikes += 1
                    neuron.firing_rate = neuron.nspikes / self.t
                    
                    neuron.potential = V[i]
                                    
                    if neuron.firing_rate < self.tfr:
                        neuron.size *= (1 + self.growth_coeff)
                    elif neuron.firing_rate >= self.tfr: 
                        neuron.size *= (1 - self.growth_coeff)
                    if neuron.size > self.size_range[1]: neuron.size = self.size_range[1]
                    if neuron.size < self.size_range[0]: neuron.size = self.size_range[0]
                

                self.update_connectivity_matrix()
                
                # Updating global lists (keep track of all attributes for every iteration)
                
                firing_rates = [neuron.firing_rate for neuron in self.neurons]
                neuron_sizes = [neuron.size for neuron in self.neurons]
                self.mfr = np.mean(np.array(firing_rates))
                
                if mfr:
                    self.tfr = self.mfr
                    
                V_list.append(V.copy())
                spikes_list.append(spikes.copy())
                firing_rates_list.append(firing_rates.copy())
                size_list.append(neuron_sizes.copy())                         
            
        print('mean firing rate:', self.mfr)    
        return spikes_list, V_list, firing_rates_list, size_list
        
        
    def compute_avalanche_sizes(self, spikes_list, avalanche_threshold, tau_seconds):

        num_steps = len(spikes_list)

        tau_steps = int(tau_seconds / 0.1)

        avalanche_sizes = []
        current_avalanche_size = 0

        for t in range(num_steps):

            window_start = max(0, t - tau_steps)
            window_end = min(num_steps, t + 1)

            n_neurons_spiking = np.sum(spikes_list[window_start:window_end])

            if n_neurons_spiking >= avalanche_threshold:
                current_avalanche_size += n_neurons_spiking
                
            else:
                if t == num_steps:
                    avalanche_sizes.append(current_avalanche_size)
                    current_avalanche_size = 0
                    
                
                elif current_avalanche_size > 0:
                    avalanche_sizes.append(current_avalanche_size)
                    current_avalanche_size = 0

        return avalanche_sizes
    
class Graphics:
    
    def draw(self, network):
        fig, ax = plt.subplots(figsize=(10, 10))
        for neuron in network.neurons:
            circle = plt.Circle(neuron.position, neuron.size, color=neuron.color, alpha=0.5)
            ax.add_artist(circle)
        ax.set_xlim(0, network.sandbox_size[0])
        ax.set_ylim(0, network.sandbox_size[1])
        ax.set_aspect('equal')
        ax.axis('off')

    def plot_neuron_activity(self, spikes_list, V_list):
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.imshow(1-np.array(spikes_list).T, aspect='auto', cmap = 'binary')
        plt.xlabel('Time')
        plt.ylabel('Neuron Index')
        plt.title('Neuron Firings Over Time')
        plt.subplot(1, 2, 2)
        plt.imshow(np.array(V_list).T, aspect='auto', cmap='winter')
        plt.xlabel('Time')
        plt.ylabel('Neuron index')
        plt.title('Membrane Potential Over Time')
        plt.colorbar(label='Membrane Potential')
        plt.tight_layout()
        plt.show()
        
    def plot_M(self, network):
        plt.figure()
        plt.imshow(network.M)
        colorbar = plt.colorbar()
        colorbar.set_label('Synaptic Connection Strength')
        plt.xlabel('Neuron index')
        plt.ylabel('Neuron index')
        plt.title('Connectivity Matrix M')

    def plot_firing_rates(self, network, firing_rates_list):
        plt.figure()
        fr_img = plt.imshow(np.array(firing_rates_list).T, aspect='auto')
        plt.colorbar(fr_img)
        fr_img.set_clim(0, 1.5 * network.tfr)
        plt.title(f'Neuron firing rates, tfr = {network.tfr}')
        plt.xlabel('Iteration #')
        plt.ylabel('Neuron index')
        return fr_img
        
    def plot_neuron_sizes(self, network, size_list):
        plt.figure()
        sz_img = plt.imshow(np.array(size_list).T, aspect='auto')
        plt.colorbar(sz_img)
        sz_img.set_clim(0, network.size_range[1])
        plt.title(f'Neuron radii, size range: {network.size_range}')
        plt.xlabel('Iteration #')
        plt.ylabel('Neuron index')
        
    def bin_avalance_durations(self, avalanche_sizes):
        if not avalanche_sizes:
            print('No avalanches detected; skipping avalanche histogram.')
            return
        plt.hist(avalanche_sizes, bins=range(1, max(avalanche_sizes) + 1), align='left', edgecolor='black')
        plt.xlabel('Avalanche Size')
        plt.ylabel('Frequency')
        plt.title('Histogram of Avalanche Sizes')
        plt.show()
        
        
    def plot_avalanche_sizes(self, network, spikes_list, min_threshold, window_size):
        
        avalanche_sizes = network.compute_avalanche_sizes(spikes_list, min_threshold, window_size)
        if not avalanche_sizes:
            print('No avalanches detected for the selected threshold/window; skipping avalanche plots.')
            return

        plt.figure(figsize=(8, 6))
        plt.hist(avalanche_sizes, bins=range(1, max(avalanche_sizes) + 1), align='left', edgecolor='black')
        plt.xlabel('Avalanche Size')
        plt.ylabel('Frequency')
        plt.title('Histogram of Avalanche Sizes')
        plt.grid(True)
        plt.show()


        hist, bin_edges = np.histogram(avalanche_sizes, bins='auto', density=True)
        if len(hist) == 0 or len(bin_edges) < 2:
            print('Insufficient avalanche data for log-log plot.')
            return

        # Find the midpoint of each bin for plotting purposes
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        valid = (bin_centers > 0) & (hist > 0)
        if not np.any(valid):
            print('Avalanche histogram has no positive bins for log-log plot.')
            return

        # Plot the histogram on a log-log scale
        plt.figure(figsize=(8, 6))
        plt.loglog(bin_centers[valid], hist[valid], marker='o', linestyle='None', color='b', label='Avalanche Sizes')

        # Fit a power-law function to the data (optional)
        try:
            fit = powerlaw.fit(avalanche_sizes)
            alpha = fit[0]  # Extract the fitted exponent

            # Plot the power-law fit
            pdf_x = bin_centers[valid]
            pdf_y = powerlaw.pdf(pdf_x, alpha)
            pdf_valid = np.isfinite(pdf_y) & (pdf_y > 0)
            if np.any(pdf_valid):
                plt.loglog(pdf_x[pdf_valid], pdf_y[pdf_valid], color='r', linestyle='--', label=f'Power-law fit (alpha={alpha:.2f})')
        except Exception as e:
            print(f'Power-law fit skipped: {e}')

        plt.xlabel('Avalanche Size')
        plt.ylabel('Probability Density')
        plt.title('Log-log plot of Avalanche Sizes')
        plt.legend()
        plt.grid(True)
        plt.show()
        
        
def main():
    network1 = Network(N=100, sandbox_size = (2048,2048), tfr = 0.01, size_range=(10,50))
    graphics = Graphics()

    network1.init_neurons()
    graphics.draw(network1)

    network1.tfr = 5 #setting the target firing rate to 5 hz
    spikes_list, V_list, firing_rates_list, size_list = network1.simulate(duration = 500, dt = 0.1, spike_input=10, input_noise=0.5, leak=0.01, animate=True)
    graphics.plot_avalanche_sizes(network1, spikes_list, 1, 10)

if __name__ == "__main__":
    main()
