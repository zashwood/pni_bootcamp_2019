import numpy as np
from scipy import stats, integrate, interpolate
from scipy.signal import convolve2d
from scipy.optimize import minimize

#############################################################################################################
# COMPUTING #################################################################################################

def field_2D(x, y, mu, sigma, R):
    # placefield_2D(...) : For locations x and y (1D or 2D arrays acceptable), return the mean response
    # (i.e. the placefield) of the neuron with center mu, width sigma, and max response R. Simulated 
    # placefields are spherical Gaussians
    #
    #     x : n x m grid of x values
    #     y : n x m grid of y values
    #    mu : 2D placefield center
    # sigma : placefield width; can be scalar, 2D, or 2 x 2
    #     R : (scalar) maximum response of neuron
    
    if sigma.size == 1:
        ## Spherical Gaussian
        P = R * np.exp(-( (x-mu[0])**2 / (2*sigma**2) + (y-mu[1])**2 / (2*sigma**2) ))
    elif sigma.size == 2:
        ## Diagonal Gaussian
        P = R * np.exp(-( (x-mu[0])**2 / (2*sigma[0]**2) + (y-mu[1])**2 / (2*sigma[1]**2) ))
    return P

def KL(p, q):
    # KL(...) : Computes KL divergence (statistical distance) between two distributions p and q, which in our
    # case take values on a grid. Flatten them and return the KL divergence, which we use for its equivalence
    # to maximum likelihood estimates of the parameters estimating our parametric approximation q of p.
    #
    #     p : n x m grid of true (empirical) posterior probabilities of a spike at that location on the grid
    #     q : n x m grid of estimated posterior probabilities
    
    return stats.entropy(p.flatten(), q.flatten())

def fit_2DGaussian(sample, x, y, v0):
    # fit_2DGaussian(...) : Fits a 2D Gaussian to the sample data given by 'sample' on the grid defined jointly
    # by x and y. 
    #
    # sample : n x m grid of true (empirical) posterior probabilities of a spike at that location on the grid
    #      x : n x m grid of x values
    #      y : n x m grid of y values
    #      v : 4-D vector of initialization (guess) parameters for optimization [mean_x mean_y width_x width_y]
    
    obj = lambda v: KL(sample, field_2D(x, y, v[0:2], v[2:], 1))
    v_fit = minimize(obj, v0).x
    return field_2D(x, y, v_fit[0:2], v_fit[2:], 1)
    
def bayesianupdate(prior, kern, spikes, spike_fields):
    # locationupdate(...) : Given a grid of probabilities of being in each location, update each by the probability
    # of moving according to the random walk, then update further by any spikes observed (or not)
    #
    #        prior : n x m grid of (prior) probabilities of being in a given location
    #         kern : 3 x 3 grid of probabilities of moving in each of 8 cardinals directions + not moving
    #       spikes : N x 1, spike raster FOR THIS UPDATE'S TIME POINT ONLY ***********
    # spike_fields : ny x nx x N array, inferred (estimated) (posterior) probability of observing a spike in each location

    N = spike_fields.shape[-1]
    posterior = convolve2d(prior, kern, mode="same", boundary="symm") ## conv2d in MATLAB can't do this ;)
    ## This single line is equivalent to this nightmare below:
    # move_prior = 0.25*prior; ## Probability of no change is 0.25
    # ## Add probabilities of cardinal moves (1/8 each)
    # move_prior(1:ny-1,:) = move_prior(1:ny-1,:) + prior(2:ny,:)/8;
    # move_prior(2:ny,:) = move_prior(2:ny,:) + prior(1:ny-1,:)/8;  
    # move_prior(:,1:nx-1) = move_prior(:,1:nx-1) + prior(:,2:nx)/8;
    # move_prior(:,2:nx) = move_prior(:,2:nx) + prior(:,1:nx-1)/8;
    # ## Add probabilities of diagonal moves (1/16 each)
    # move_prior(1:ny-1,1:nx-1) = move_prior(1:ny-1,1:nx-1) + prior(2:ny,2:nx)/16;
    # move_prior(2:ny,1:nx-1,:) = move_prior(2:ny,1:nx-1,:) + prior(1:ny-1,2:nx)/16;
    # move_prior(1:ny-1,2:nx) = move_prior(1:ny-1,2:nx) + prior(2:ny,1:nx-1)/16;
    # move_prior(2:ny,2:nx) = move_prior(2:ny,2:nx) + prior(1:ny-1,1:nx-1)/16;
    # ## Add additional probabilities of not moving on edge and corner sites
    # move_prior(1,:) = move_prior(1,:) + 0.25*prior(1,:);
    # move_prior(ny,:) = move_prior(ny,:) + 0.25*prior(ny,:);
    # move_prior(:,1) = move_prior(:,1) + 0.25*prior(:,1);
    # move_prior(:,nx) = move_prior(:,nx) + 0.25*prior(:,nx);
    # ## Remove double counting for corner sites
    # move_prior(1,1) = move_prior(1,1) - prior(1,1)/16;
    # move_prior(ny,1) = move_prior(ny,1) - prior(ny,1)/16;
    # move_prior(1,nx) = move_prior(1,nx) - prior(1,nx)/16;
    # move_prior(ny,nx) = move_prior(ny,nx) - prior(ny,nx)/16;
    
    ## Now calculate updated probability based on spikes or no spikes
    ## ==> If a cell spikes, the decoded location drifts strongly towards that neuron's spike field
    for cell in range(N):
        if spikes[cell]: ## Spike
            posterior = posterior * spike_fields[:,:,cell]
        else:              ## No spike
            posterior = posterior * (1 - spike_fields[:,:,cell])

    posterior /= np.sum(posterior)
    prior = posterior      ## New prior for the next time-bin
    
    return posterior, prior

#############################################################################################################
# PLOTTING ##################################################################################################

def downsample_for_plot(data, rate=100):
    # downsample_for_plot(...) : Downsample the spike data at a rate:1 samples along the given axis. 
    # Spike rasters usually have too many sample to be visible in plots. Use this to fix that.
    #
    #   data : nt x N dataset of spiking activity (1 and 0) over nt time samples
    #   rate : (scalar) downsample rate, # samples in --> 1 sample out

    filt = np.ones((1,rate)) ## This filter counts spikes every [rate] samples
    blurdata = convolve2d(data, filt, mode="same")
    return blurdata[:,::rate]

def refreshfig(fig, posterior, point, limits):
    # refreshfig(...) : Refresh a figure in a loop with an animation-like effect. Plots a point on top of an image
    # representing the (posterior) uncertainty of that point's location.
    #
    #       fig : handle to the figure
    # posterior : ny x nx grid of inferred probability of being in position (x,y), i.e. estimate of location
    #     point : 2 x 0 array of the true position (x*,y*)
    #    limits : 4 x 0 array of x-y limits for figure, formatted specifically for the "extent" argument to plt.imshow()
    
    ## Refresh "animated" figure during the decoding computation
    plt.clf()
    ## Decoded (inferred) probability of position
    plt.imshow(posterior, cmap=plt.cm.bone, origin="lower", extent=limits) 
    ## Superpose the actual position
    plt.plot(point[1], point[0], "or")
    fig.canvas.draw()
    fig.canvas.flush_events()
    
    