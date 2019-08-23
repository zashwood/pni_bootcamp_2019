
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy import stats, integrate, interpolate
from scipy.signal import convolve2d
from scipy.optimize import minimize

import decoding_utils as utils

def groundtruth_placefields(xgrid, ygrid, place_centers, widths, rmax):
    # groundtruth_placefields(...) : Given a 2D grid of x-y coordinates and a list of place field params,
    # generate the corresponding place fields and return a plot of at most 4 of them.
    #          xgrid : n x m grid of x values
    #          ygrid : n x m grid of y values 
    #  place_centers : N x 2 list of 2D place field centers
    #         widths : N x 1 list of place field widths (ground truth assumed spherical)
    #           rmax : N x 1 list of place field maximum mean firing rates

    place_fields = np.zeros((xgrid.shape[0],xgrid.shape[1],rmax.size)) ## Empty array to fill with place-fields
    
    for cell in range(rmax.size):
        # Compute scaled spherical Gaussian, centered on the x- and
        # y-coordinates of the place-field center for each cell
        place_fields[:, :, cell] = utils.field_2D(xgrid, ygrid, place_centers[cell,:], widths[cell], rmax[cell]) 

    return place_fields 
    
def spiking_randomwalk_grid(nt, dt, nx, ny, place_fields, step_prob=[0.25, 0.5, 0.25]):
    # spiking_randomwalk_grid(...) : Simulate a 2D random walk with variable probabilities of motion on an integer
    # grid of size nx x ny for nt time steps. Generate spikes from location according to neurons in place_fields.
    #             nt : number of time steps
    #             dt : time increment per time step, for spike rates
    #             nx : size of xgrid (0 -- nx)
    #             ny : size of ygrid (0 -- ny)
    #   place_fields : ny x nx x N set of placefields for each neuron
    
    N = place_fields.shape[-1]
    spikes = np.zeros((N,nt)) ## Record presence or absence of spike as 1 or 0
    position = np.zeros((nt,2))                    ## To store x- and y-coordinates of actual position
    position[0] = [nx/2,ny/2]                      ## Initial position in center (could generalize this...)
    
    step_set = [-1,0,1]
    for i in range(1,nt):
        step = np.random.choice(a=step_set, size=2, p = step_prob) 

        ## Enforce boundary conditions
        newx = np.maximum(np.minimum(position[i-1,1]+step[1],nx-1),0)
        newy = np.maximum(np.minimum(position[i-1,0]+step[0],ny-1),0)

        position[i] = [newy, newx]

        ## Produce a spike with probability of rate*dt for each neuron
        spikes[:,i] = (np.random.rand(N) < dt*np.squeeze(place_fields[int(newy),int(newx),:]))
    return position, spikes
    
def estimate_spikefields(position, spikes, xgrid, ygrid):
    # estimate_spikefields(...) : Compute empirical posterior probability of a spike for each location. Do this
    # by measuring the proportion of (# spikes in each location) to (# timesteps spent in each location).
    # Use a parametric fit to fit a 2D isotropic Gaussian to this empirical posterior to estimate the spikefields.
    #       position : nt x 2 integer array of grid positions over time
    #         spikes : N x nt spike raster
    #          xgrid : ny x nx x-coordinate grid
    #          ygrid : ny x nx y-coordinate grid
    
    nt, nx, ny, N = spikes.shape[1], xgrid.shape[1], xgrid.shape[0], spikes.shape[0]
    
    ## Count time spent at each position
    thalf  = int(np.round(nt/2))   ## Time-point ending first half of data
    xybase = np.zeros((ny,nx))  ## To record time spent at a given (x,y) coordinate
    for i in range(thalf):
        xi = position[i,1]
        yi = position[i,0]
        xybase[yi,xi] += 1

    ## Count spikes produced for each "neuron" at each position
    xyhist = np.zeros((ny,nx,N))
    spike_fields = np.zeros((ny,nx,N))
    prob_spike = np.zeros(N)

    for cell in range(N):
        ispk = np.where(spikes[cell,0:thalf])[0] ## Extract the time-points of spikes

        for i in range(ispk.size):
            xi = position[ispk[i],1] ## x-coordinate at spike-time
            yi = position[ispk[i],0] ## y-coordinate at spike-time
            xyhist[yi,xi,cell] += 1

        ## Divide by time spent at each location, to get probability of a spike in each time bin when at that location.
        ## \propto Pr(spike | location)
        xyhist[:,:,cell] /= np.maximum(xybase,1)

        ## Initialize parameters of Gaussian fit
        mu0 = np.array([np.sum(xgrid * xyhist[:,:,cell]) / np.sum(xyhist[:,:,cell]), \
                        np.sum(ygrid * xyhist[:,:,cell]) / np.sum(xyhist[:,:,cell])])
        sigma0 = np.array([10,10])
        ## Fit
        spike_fields[:,:,cell] = utils.fit_2DGaussian(xyhist[:,:,cell], xgrid, ygrid, np.concatenate((mu0,sigma0)))

        ## Measure average prob of spike for posterior computation downstream, Pr(spike)
        prob_spike[cell] = np.mean(xyhist[:,:,cell])

    ## Assume prob of location is uniform, Pr(location)
    prob_loc = 1/(nx+ny)
    
    return spike_fields, prob_spike, prob_loc, xyhist

def estimate_placefields(position, spikes, dt, xgrid, ygrid):
    # estimate_placefields(...) : Compute and renormalize the spike fields to estimate the place fields.
    #       position : nt x 2 integer array of grid positions over time
    #         spikes : N x nt spike raster
    #          xgrid : ny x nx x-coordinate grid
    #          ygrid : ny x nx y-coordinate grid
    
    nt, nx, ny, N = spikes.shape[1], xgrid.shape[1], xgrid.shape[0], spikes.shape[0]
    spike_fields, prob_spike, prob_loc, xyhist = estimate_spikefields(position, spikes, xgrid, ygrid)
    
    nospike_fields = np.zeros(spike_fields.shape)
    est_place_fields = np.zeros(spike_fields.shape)

    for cell in range(N):

        ## spike_fields is the probability of a spike when in a particular location, so should be normalized 
        ## -- here by multiplying by (nx*ny) we assume implicitly that the probability of being at any location is 
        ## equal at 1/(nx*ny) so that we have used Bayes' Theorem:
        ## P(spike | location) = P(location | spike) * P(spike) / P(location)
        spike_fields[:,:,cell] *= (prob_spike[cell]*(nx*ny)/np.sum(spike_fields[:,:,cell]))

        ## Convert spike probability per time bin to a rate for estimated place-field
        est_place_fields[:,:,cell] = spike_fields[:,:,cell]/dt
        
    return spike_fields, est_place_fields, xyhist

def decode(xx, yy, position, spikes, spike_fields, step_prob=[0.25, 0.5, 0.25], plot_flag=False):
    # decode(...) : Use Bayesian updating to translate estimated place fields to update estimates 
    # of positions given spikes and given current position estimates.
    #             xx : nx x 1 array of x-locations
    #             yy : ny x 1 array of y-locations
    #       position : nt x 2 array of locations in the space over time
    #         spikes : N x nt, spike raster
    #   spike_fields : ny x nx x N array, inferred (estimated) (posterior) probability of observing a spike in each location
    #      step_prob : 3 x 1 probability of moving (assumed symmetric, could be more general if needed...)
    #      plot_flag : bool, if True produces an animated plot of the decoder
    
    nt, nx, ny, N = position.shape[0], spike_fields.shape[1], spike_fields.shape[0], spike_fields.shape[2]
    thalf = int(nt/2)

    prior = np.ones((ny,nx))/(nx*ny) ## Uniform prior
    ## We will use the previous probability estimate of position to produce a prior on the 
    ## current probability distribution, taking into account the (known) probability of moving
    ## in different directions.
    ## On the grid, there are 9 possible new positions given the old position:
    ##  UL (1/16) U (1/8)  UR(1/16)
    ##   L (1/8)  0 (1/4)  R (1/8)
    ##  DL (1/16) D (1/8)  DR(1/16)
    kern = np.expand_dims(step_prob, axis=0)*np.expand_dims(step_prob, axis=1)

    move_prior = np.zeros((ny,nx))       ## To be updated each time-step
    xyest = np.zeros((nt-thalf, 2))      ## Estimates of x- and y-coordinates

    fig, ax = plt.subplots()
    for i in range(thalf,nt):

        posterior, prior = utils.bayesianupdate(prior, kern, spikes[:,i], spike_fields)

        ## Animated plot that updates every 100 steps if flagged
        if plot_flag:
            if np.mod(i,100) == 0:
                utils.refreshfig(fig, posterior, position[i,:], [xx[0],xx[-1],yy[0],yy[-1]])

        ## Decode a point estimate of the x-coordinate and y-coordinate for later performance statistics 
        xyest[i-thalf,0] = xx[np.argmax(np.sum(posterior, axis=1))]
        xyest[i-thalf,1] = yy[np.argmax(np.sum(posterior, axis=0))]
        
    return xyest

def compute_scores(xyest, position):
    # decode(...) : Calculate a couple of scores to assess how well position was tracked 
    #          xyest : nt x 2 array of decoded position over time
    #       position : nt x 2 array of actual position over time
    
    ## Correlation
    yxc = np.corrcoef(xyest[:,0], position[:,0])[1,0]
    xxc = np.corrcoef(xyest[:,1], position[:,1])[1,0]
    score = (yxc+xxc)/2

    ## Mean-squared error
    MSE = np.mean( np.sum((xyest - position)**2, axis=1) )
    
    return score, MSE

## -------------------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------------------

if __name__ == '__main__': ## This code only executes if you call this script directly from the terminal/command prompt.
    print("0")