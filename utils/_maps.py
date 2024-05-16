import numpy as np
import astropy as ap
import astropy.units as u
from astropy.coordinates import SkyCoord, Galactic

from astropy_healpix import HEALPix
import healpy as hp

from matplotlib.path import Path

from scipy.spatial import Delaunay

# ======================================================================
# ======================================================================
                    # MAP GENERATION FUNCTIONS
    # References:
        # Development Directory: od/projects/pswavelets/gce/matt_post_energy
        # gce/matt_post_energy/skymap_energy.ipynb
# ======================================================================
# ======================================================================

def generate_energy_bins_(choice=1):
    """
    Generate MeV energy bins 
        Settings determined by Matt
        
    :param choice = choice for energy bins
    
    :output = list containing
        [0] = list containing lower bound of energy bins
        [1] = list containing the centers of each energy bin
    """
    if choice == 0:
        E_min = 1e2 #min gamma-ray energy (MeV)
        E_max = 1e9 #max gamma-ray energy (MeV)
        E_factor = 1.77828 # gamma-ray energy factor

        E = E_min
        energy_list = []

        while E <= E_max:
            energy_list.append(E)
            E *= E_factor
        energy_list.append(E)
        energy_centers = np.concatenate([energy_list,[energy_list[-1]*1.7]])
        energy_centers = 0.5*(energy_centers[0:-1]+energy_centers[1:])
        return [energy_list, energy_centers]
    
    elif choice == 1:
        #this is the lower and upper bounds of the energy bins
        energy_bins = np.array([43.8587, 57.0013,74.082,96.2812,125.133,162.629,211.362,274.698, 
                                357.014,463.995,603.034,783.737,1018.59,1323.82,1720.51,2236.07,
                                2906.12,3776.96,4908.75,6379.69,8291.4,10776.0,14005.1,18201.8,
                                23656.1,30744.8,39957.6,51931.2,67492.7,87717.4,114002.0,148164.0,
                                192562.0,250265.0,325258.0,422724.0,549396.0,714027.0,927989.0])
        energy_list = energy_bins[:-1]
        
        #bin centers
        energy_centers = np.logspace( np.log10(50),np.log10(814008.272176),38)
        return [energy_list, energy_centers]

def psf_(energy,popt): 
    """
    Gaussian fit of PSF:
        x = log(E)
        PSF(E) = ax^2 + bx + c
    
    :param energy = energy (MeV)
    :param popt = quadratic fit parameters [a,b,c] 
    
    :output = psf value (rad)
    """
    x = np.log(energy)
    poly = np.poly1d(popt)
    y = poly(x)
    p = np.exp(y) #degrees
    return p*np.pi/180 #radians

# map generation
def healpix_counts_to_events_(random_map, energy_list, popt, nside=256):
    """
    Convert a map of counts on healpix pixels to events.
    
    :param random_map = random generation of (background) map
    :param energy_list = list specifying energy bins (see generate_energy_bins_ fct)
    :param popt = quadratic fit parameters (see psf_ fct)
    :param nside = number of sides of healpix pixels (default = 256)
    
    :output = dictionary with event data
        :key 'coords' = true angular coordinates [l,b] of events (rad)
        :key 'energies' = energies of events (MeV)
        :key 'smeared_coords' = smeared angular coordinates [l,b] of events (rad)
    """
    #theta and phi locations of each cell
    theta_list,phi_list = hp.pix2ang(nside, range(hp.nside2npix(nside)))
    res = hp.nside2resol(nside) 
    
    randomized_dict = {}
    # generate unsmeared events
    for ie in range(len(energy_list)):
        energy_dict = {}

        theta_events = np.array([])
        phi_events = np.array([])

        low = energy_list[ie]
        if ie!= len(energy_list)-1:
            high = energy_list[ie+1]
        else:
            #high = 1.8*energy_list[-1]          model: SA0
            high = 927989.0                     # model: ilias_60x60

        for i in range(len(random_map[ie])):
            counts = np.random.poisson(lam=random_map[ie,i])
            theta = theta_list[i]
            phi = phi_list[i]

            #generate 'counts' random numbers centered on theta,phi, smeared by a gaussian
            # NEED TO MODIFY THIS (SEE SLACK 04/27)
            if counts==0:
                continue
            else:
                theta_events = np.concatenate([theta_events,theta+np.random.normal(size=int(counts),scale=res)])
                phi_events = np.concatenate([phi_events,phi+np.random.normal(size=int(counts),scale=res)])

        #correctly orient in Galactic coords
        phi_events = -phi_events
        theta_events = theta_events-np.pi/2 
        theta_events = -theta_events 
        coordinates = SkyCoord(l=phi_events*u.rad,b=theta_events*u.rad,frame=Galactic)        
        
        # generate energy events
        energy_events = np.random.uniform(size=len(coordinates),low=low,high=high)       
        
        # smear events with psf
        theta_events = coordinates.b.rad
        phi_events = coordinates.l.rad

        angular_offsets = np.random.normal(size=len(theta_events),scale=psf_(energy_events, popt))*u.rad
        azimuthal_offsets = np.random.uniform(size=len(theta_events),low=0,high=2*np.pi)*u.rad

        coordinates_new = coordinates.directional_offset_by(azimuthal_offsets, angular_offsets)
        
        # update dictionaries
        energy_dict['coords'] = coordinates
        energy_dict['energies'] = energy_events
        energy_dict['smeared_coords'] = coordinates_new
        randomized_dict[ie] = energy_dict
    return randomized_dict

def beta_healpix_counts_to_events_(random_map, energy_list, popt, nside=256):
    """
    (Beta)
    Same as the above. However, the sampling is slightly more accurate. The above samples over phi 
    over the same scales for all theta. This means the scale really depends on theta. 
    
    To do what we intend, sample from a Gaussian over the pixels, we use the same technique to 
    smear coordinates with the psf.
    """
    #theta and phi locations of each cell
    theta_list,phi_list = hp.pix2ang(nside, range(hp.nside2npix(nside)))
    
    #correctly orient in Galactic coords
    phi_list = -phi_list
    theta_list = theta_list-np.pi/2 
    theta_list = -theta_list 

    # size of pixels
    res = hp.nside2resol(nside) 
    
    randomized_dict = {}
    # generate unsmeared events
    for ie in range(len(energy_list)):
        energy_dict = {}

        theta_events = np.array([])
        phi_events = np.array([])

        low = energy_list[ie]
        if ie!= len(energy_list)-1:
            high = energy_list[ie+1]
        else:
            high = 1.8*energy_list[-1]

        for i in range(len(random_map[ie])):
            counts = np.random.poisson(lam=random_map[ie,i])
            theta = theta_list[i]
            phi = phi_list[i]

            #generate 'counts' random numbers centered on theta,phi, smeared by a gaussian
            if counts==0:
                continue
            else:
                # pixel coordinate as a SkyCoord object
                pix_coordinate = SkyCoord(l=phi*u.rad,b=theta*u.rad,frame=Galactic)

                #generate offsets from point source with gaussian with appropriate psf 
                angular_offsets = np.random.normal(size=int(counts),scale=res)*u.rad
                azimuthal_offsets = np.random.normal(size=int(counts),scale=res)*u.rad

                # sample events as offset from pixel coordinate
                coordinates = pix_coordinate.directional_offset_by(azimuthal_offsets, angular_offsets)
                phi_events = np.concatenate((phi_events, coordinates.l.rad)) # combine as np array
                theta_events = np.concatenate((theta_events, coordinates.b.rad)) # combine as np array
            
        coordinates = SkyCoord(l=phi_events*u.rad,b=theta_events*u.rad,frame=Galactic) # write as SkyCoord obj
        
        # generate energy events
        energy_events = np.random.uniform(size=len(coordinates),low=low,high=high)       
        
        # smear events with psf
        theta_events = coordinates.b.rad
        phi_events = coordinates.l.rad

        angular_offsets = np.random.normal(size=len(theta_events),scale=psf_(energy_events, popt))*u.rad
        azimuthal_offsets = np.random.uniform(size=len(theta_events),low=0,high=2*np.pi)*u.rad

        coordinates_new = coordinates.directional_offset_by(azimuthal_offsets, angular_offsets)
        
        # update dictionaries
        energy_dict['coords'] = coordinates
        energy_dict['energies'] = energy_events
        energy_dict['smeared_coords'] = coordinates_new
        randomized_dict[ie] = energy_dict
    return randomized_dict

def generate_pointsource_(l,b,photon_number,energy_list,popt,energy_bin=0,randomize_number=False):
    """
    Generate point source events from central position and Poisson average.
        * The total number of events are drawn from a Poisson distribution at each energy bin.
        * Distribution of events is uniform in energy.
    
    :param l = longitude of point source (l \in [0,2\pi])
    :param b = latitude of point source (b \in [-\pi/2, \pi/2])
    :param energy_list = list specifying energy bins (see generate_energy_bins_ fct)
    :param popt = quadratic fit parameters (see psf_ fct)
    :param energy_bin = index corresponding to energy bin of choice
    :param randomize_number = boolean type indicating whether number of events drawn randomly
    
    :output = dictionary with point source data
        :key 'center_coords' = (l,b)-coordinates of point source center (rad) 
        :key 'energies' = energies of point source events (MeV)
        :key 'coords' = true (l,b)-coordinates of point source events (rad) 
        :key 'smeared_coords' = true (l,b)-coordinates of point source events (rad) 
    """
    #l, b in rad
    if randomize_number:
        counts = np.random.poisson(lam=photon_number)
    else:
        counts = photon_number
    
    #generate energies appropriate for the energy bin (uniform distribution for now)
    low = energy_list[energy_bin]
    if energy_bin!= len(energy_list)-1:
        high = energy_list[energy_bin+1]
    else:
        high = 1.8*energy_list[-1]

    energy_events = np.random.uniform(size=counts,low=low,high=high)    
    
    #create list of skycoords at the ps center
    l_list = l*np.ones(counts)
    b_list = b*np.ones(counts)
    
    coordinates = SkyCoord(l=l_list*u.rad,b=b_list*u.rad,frame=Galactic)
    
    #generate offsets from point source with gaussian with appropriate psf 
    angular_offsets = np.random.normal(size=counts,scale=psf_(energy_events, popt))*u.rad
    azimuthal_offsets = np.random.uniform(size=counts,low=0,high=2*np.pi)*u.rad
    
    
    coordinates_new = coordinates.directional_offset_by(azimuthal_offsets, angular_offsets)
    
    ps_dict = {}
    ps_dict['center_coords'] = np.array([l,b]) 
    ps_dict['energies'] = energy_events
    ps_dict['coords'] = coordinates
    ps_dict['smeared_coords'] = coordinates_new
    
    return ps_dict

def generate_point_sources_inside_curve_(N_ps, N_counts, npix, lon_edge, lat_edge, energy_list, popt,energy_bin=0,randomize_number=True):
    """
    Sample point sources uniformly on a curve. (Development: random_sampling_spherical_section.ipynb)
        * The total number of events are drawn from a Poisson distribution at each energy bin.
        * Distribution of events is uniform in energy.
    
    :param N_ps = number of point sources
    :param N_counts = average number of counts per point source
    :param npix = pixel number
    :param lon_edge = longitude of father pixel edge
    :param lat_edge = latitude of father pixel edge
    :param energy_list = list specifying energy bins (see generate_energy_bins_ fct)
    :param popt = quadratic fit parameters (see psf_ fct)
    :param energy_bin = index corresponding to energy bin of choice
    :param randomize_number = boolean type indicating whether number of events drawn randomly
    
    :output = list of dictionary with point source data
        dict structure: [index][key]
            :index = point source label
            :key 'center_coords' = (l,b)-coordinates of point source center (rad) 
            :key 'energies' = energies of point source events (MeV)
            :key 'coords' = true (l,b)-coordinates of point source events (rad) 
            :key 'smeared_coords' = true (l,b)-coordinates of point source events (rad) 
        
    NOTE: Does not work for disconnected pixels. Need to carefully derive that result separately if necessary.
          Alternative choice is to sample them in coordinates where the pixel is not disconnected.
          I.e., shift by l by +np.pi. Sample. Then, shift back.
    """  
    
    # define "rectangular" region containing the pixel
    l_edge = lon_edge - np.pi ; b_edge = lat_edge 

    l_edge_min = np.min(l_edge) ; l_edge_max = np.max(l_edge) 
    b_edge_min = np.min(b_edge) ; b_edge_max = np.max(b_edge)

    delta_l = l_edge_max - l_edge_min
    sin_b_min = np.sin(b_edge_max)
    delta_sin_b = - np.sin(b_edge_max) + np.sin(b_edge_min)

    ps_dict_list = []
    l_list = []
    b_list = []

    # sample points in rectangular region containing pixel and keep ones inside pixel
    while len(ps_dict_list) < N_ps:
        l = delta_l * np.random.random(1) + l_edge_min
        b = ( np.pi/2 - np.arccos( (delta_sin_b * np.random.random( 1 ) + sin_b_min )) ) 

        ps_dict = generate_pointsource_(l,b,N_counts,energy_list,popt,energy_bin,randomize_number) 
        l_ps = ps_dict['smeared_coords'].l.rad
        b_ps = ps_dict['smeared_coords'].b.rad

        # after loading the data, we need our angular coordinates to be given by
        # longitude ([0,2\pi]) and latitude (-\pi, \pi)
        phi_ps = l_ps.copy()
        phi_ps[phi_ps>np.pi] = phi_ps[phi_ps>np.pi]-2*np.pi

        lon_ps = phi_ps + np.pi
        lat_ps = b_ps

        lon_within, lat_within = find_points_inside_pixel_(lon_ps, lat_ps, npix, lon_edge, lat_edge) # 2D array
        if len(lon_within) == len(lon_ps):
            ps_dict_list.append(ps_dict)
            l_list.append(l)
            b_list.append(b)
        else:
            continue
        
    ps_loc = np.stack((l_list, b_list), axis = -1)
    ps_loc = np.squeeze(ps_loc, axis = 1)
    
    return [ps_dict_list, ps_loc]

def dict_to_array_(randomized_dict,energy_list):
    """
    Convert map data from dictionary to array.
    
    :param randomized_dict = dictionary containing data of randomized map (see healpix_counts_to_events_ fct)
    :param energy_list = list specifying energy bins (see generate_energy_bins_ fct)
    
    :output = Array of datapoints
        axis 0 = index of datapoint
        axis 1 = properties of datapoint
            [0] = longitude (rad)
            [1] = latitude (rad)
            [2] = energy (MeV)
            [3] = energy bin
    """
    # produce array of datapoints 
    ## index 0: event number
    ## index 1: event properties
    ### 0: galactic longitude (rad)
    ### 1: galactic latitude (rad)
    ### 2: energy (MeV)
    ### 3: energy bin 
    all_events = []
    for ie in range(len(energy_list)):
        coords = randomized_dict[ie]['smeared_coords']
        energies = randomized_dict[ie]['energies']

        l_events = coords.l.rad
        b_events = coords.b.rad
        e_events = energies
        e_bin_type = ie * np.ones(len(l_events))

        if ie == 0:
            events = np.stack((l_events, b_events, e_events, e_bin_type), axis = -1)
        else:
            events_ie = np.stack((l_events, b_events, e_events, e_bin_type), axis = -1)
            events = np.concatenate((events, events_ie), axis = 0)
    return events

#recenter coordinates to a new coordinate system centered on a point
def recenter_(coordinates,coordinate_center):
    coordinates_new = coordinates.transform_to(SkyOffsetFrame(origin=coordinate_center)) 
    return coordinates_new
# ======================================================================
# ======================================================================


# ======================================================================
# ======================================================================
                # TANGENT PLANE PROJECTION FUNCTIONS
    # References:
        # Development Directory: od/projects/pswavelets/gce/notebooks
        # Basics: https://en.wikipedia.org/wiki/Spherical_coordinate_system
        # Projection Method: tangent_plane_projection_exp.pdf 
        # Mathematica Code: tangent_plane_projection_final.nb
# ======================================================================
# ======================================================================

def r_(r,theta, phi):
    """
    Create (x,y,z) position vector from spherical coordinates.
    
    :param r = radial distance
    :param theta = polar angle (rad)
    :param phi = azimuthal angle (rad)
    
    :output = array [x,y,z]
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.array([x,y,z])

def theta_hat_(theta,phi):
    """
    Create theta_hat unit-vector from (theta,phi)
    
    :param theta = polar angle (rad)
    :param phi = azimuthal angle (rad)
    
    :output = array [theta_hat_x, theta_hat_y, theta_hat_z]
    """
    theta_hat_x = np.cos(theta) * np.cos(phi)
    theta_hat_y = np.cos(theta) * np.sin(phi)
    theta_hat_z = - np.sin(theta)
    return np.array([theta_hat_x, theta_hat_y, theta_hat_z])

def phi_hat_(theta,phi):
    """
    Create phi_hat unit-vector from (theta,phi)
    
    :param theta = polar angle (rad)
    :param phi = azimuthal angle (rad)
    
    :output = array [phi_hat_x, phi_hat_y, phi_hat_z]
    """
    phi_hat_x = -np.sin(phi)
    phi_hat_y = np.cos(phi)
    phi_hat_z = 0
    return np.array([phi_hat_x, phi_hat_y, phi_hat_z])

def tangent_plane_proj_(lat, lon, lat_c, lon_c):
    """
    Project points on unit sphere into the plane tangent to some fixed "center" point.
        (see Projection Method reference for more detail)
    
    :param lat = latitude of point to be projected (rad)
    :param lon = longitude of point to be projected (rad)
    :param lat_c = latitude of "center" defining tangent plane (rad)
    :param lon_c = longitude of "center" defining tangent plane (rad)
    
    :output = (x,y) coordinates of point projected into the sphere
        x = longitude relative to (lon_c, lat_c) [East = Positive] (rad)
        y = latitude relative to (lon_c, lat_c) [North = Positive] (rad)
            (see snewflip2d function in Mathematica reference above)
        buf_output.shape = (N_data, 2)
        
        Note: Transpose at the end is necessary for broadcasted output to have dimension (N_data, 2)
    """
    # convert longitude [range: (0,2\pi)] to azimuthal angle [range: (-\pi,\pi)]
    # convert latitude [range: (-\pi/2, \pi/2)] to polar angle [range: (\pi,0)]
    phi = lon - np.pi 
    theta = np.pi / 2  - lat
    
    phi_c = lon_c - np.pi 
    theta_c = np.pi / 2  - lat_c
    
    return np.array([np.sin(theta) * np.sin(phi - phi_c), 
                    -np.cos(theta_c) * np.cos(phi - phi_c) * np.sin(theta)
                     + np.cos(theta) * np.sin(theta_c)]).T 

def inv_tangent_plane_proj_(x,y,lon_c,lat_c):
    """
    Project points from plane tangent to some fixed "center" point to the unit sphere.
        (see Projection Method reference for more detail)
    
    :param x = longitude relative to (lon_c, lat_c) [East = Positive] (rad)
    :param y = latitude relative to (lon_c, lat_c) [North = Positive] (rad)
    :param lat_c = latitude of "center" defining tangent plane (rad)
    :param lon_c = longitude of "center" defining tangent plane (rad)
    
    :output = (l,b)-coordinates of point projected into the sphere
        lon = longitude of point to be projected (rad)
        lat = latitude of point to be projected (rad)
            (see snewflip2d function in Mathematica reference above)
        buf_output.shape = (N_data, 2)
        
        Note: Transpose at the end is necessary for broadcasted output to have dimension (N_data, 2)
    """
    phi_c = lon_c - np.pi 
    theta_c = np.pi / 2.  - lat_c
    
    r_c = r_(1., theta_c, phi_c)
    phi_hat = phi_hat_(theta_c, phi_c)
    theta_hat = theta_hat_(theta_c, phi_c)
    
    # (x,y,z) coordinates of point on plane
    r = r_c + x * phi_hat - y * theta_hat
    # scaling factor required to translate point r down to the sphere (see reference)
    b = ( np.dot(r, r_c) - np.sqrt(np.dot(r, r_c)**2. - (np.dot(r, r) - 1.) ) )
    
    r_sph = r - b * r_c
    x_sph, y_sph, z_sph = r_sph
    theta, phi = [np.arccos(z_sph), np.arctan2(y_sph, x_sph)]
    
    lon = phi + np.pi
    lat = np.pi/2. - theta
    return np.array([lon,lat]).T

# ======================================================================
# ======================================================================

# ======================================================================
# ======================================================================
                # GROUPING POINTS INTO HEALPIX PIXELS
    # References:
        # Development Directory: od/projects/pswavelets/gce/notebooks
            # sections_of_sphere_to_plane.ipynb
            # patch_generator.ipynb
            # boundary_of_pointset.ipynb
        # Basics: The HEALPix Primer (for numbering scheme)
# ======================================================================
# ======================================================================

def healpix_edge_generator_(NSIDE=4, step=100):
    """
    Generate HEALPIX edges (easy)
        * Remove jump discontinuities from the edges associated with
          periocity of the angular coordinates (hard, done by eye)
        * Close the loops in edges
        
        Reason: Group points in (l,b) without having to project
        
    :param NSIDE = number of sides of healpix pixels (default = 4)
    :param step = determines points needed to build pixel edges (default = 100)
    
    :output = array of healpix edges
        axis 0 = pixel number
        axis 1 = index of point in edge
        axis 2 = [lon,lat]
        
        output.shape = (Npix, Ndata, 2)
    
    WARNING: Does not work for NSIDE != 4
    """
    if int(NSIDE) != NSIDE or int(step) != step: # error-check
        return print('Invalid Input')
    
    NPIX = hp.nside2npix(NSIDE) # number of father pixels
    N_edge_points = NSIDE * step + 1 # number of points in each pixel edge (extra point so loop is closed)
    arr_edge_points = np.zeros((NPIX, N_edge_points, 2)) # preallocation of edge points array

    for npix in range(NPIX):
        ahp = HEALPix(nside = NSIDE, order = 'ring')
        lon, lat = ahp.boundaries_lonlat([npix], step = step)

        lon = lon[0]
        lat = lat[0]
            
        # NSIDE = 4 is optimized so all father pixel edges are closed curves in (lon,lat)-space
        # this allows us to identify points within a father pixel before projection
        # otherwise, we need to project all points, project the edges, and identify the points then
        if NSIDE == 4: 
            lon = lon  / u.rad
            lat = np.pi/2 * u.rad - lat # convert lat to polar angle theta
                
            delta_lon = np.diff(lon)
            n_lon = np.where(delta_lon < -6)[0] + 1
            lon[n_lon] += 2 * np.pi

            n_0 = np.where(lon == 0.)[0]
            for k in n_0:
                if (lon[k] == 0.) and (np.abs(lon[k - 1] - np.pi) < 0.001):
                    while lon[k] == 0.:
                        lon[k] += np.pi
                        k += 1
                        if k == len(lon):
                            break
                elif (lon[k] == 0.) and (np.abs(lon[k + 1] - np.pi) < 0.001):
                    while lon[k] == 0.:
                        lon[k] += np.pi
                        k -= 1
                        if k == len(lon):
                            break
                elif (lon[k] == 0.) and (np.abs(lon[k - 1] - 2 * np.pi) < 0.001):
                    while lon[k] == 0.:
                        lon[k] += 2 * np.pi
                        k += 1
                        if k == len(lon):
                            break
                elif (lon[k] == 0.) and (np.abs(lon[k + 1] - 2 * np.pi) < 0.001):
                    while lon[k] == 0.:
                        lon[k] += 2 * np.pi
                        k -= 1
                        if k == len(lon):
                            break

            if lon[0] == np.pi and (np.abs(lon[1] - np.pi) > 0.1):
                lon[0] = lon[1]
            if lon[0] == 0. and (np.abs(lon[1] - 2 * np.pi) < 0.1):
                lon[0] = 2 * np.pi

            n_pi = np.where(lon == np.pi)[0]
            if npix == 191:
                if np.abs(lon[n_pi[-1] + 1] - np.pi) > 0.01:
                    k = n_pi[-1]
                    while lon[k] == np.pi:
                        lon[k] = lon[n_pi[-1] + 1]
                        k -= 1
                # extend boundary to close the sphere
                n_fix_by_eye = np.where(lon > 6.2676328)[0]
                lon[ n_fix_by_eye[:-1] ] = 2 * np.pi

            lon  = lon  * u.rad
            lat = np.pi / 2 * u.rad - lat # convert polar angle theta back to latitude

        # close loop 
        lon = np.append(lon, lon[0])
        lat = np.append(lat, lat[0])

        arr_edge_points[npix, :, 0] = lon.value
        arr_edge_points[npix, :, 1] = lat.value

    return arr_edge_points

def father_pixel_center_generator_(arr_edge_points):
    """
    Generate center of father pixels from edge arrays. 
        Center: Midpoint between max/min of l and b
        
    :param arr_edge_points = array of healpix edges (see healpix_edge_generator_ fct)
    
    :output = array of (l,b) coordinates (rad)
        output.shape = (Npix, 2)
    """
    NPIX, N_edge_points = arr_edge_points.shape[:-1] # number of father pixels, number of points in each pixel edge 
    Nhalf = int(N_edge_points / 2)
    disconnected_pixels = [40,72,104,136] 
    
    arr_c = np.zeros((NPIX,2))
    for npix in range(NPIX):
        lon = arr_edge_points[npix, :, 0]
        lat = arr_edge_points[npix, :, 1]

        if npix not in disconnected_pixels:
            lon_c = 0.5 * (np.max(lon) + np.min(lon))
            lat_c = 0.5 * (np.max(lat) + np.min(lat))
        else:
            lon_c = 0
            lat_c = 0.5 * (np.max(lat[:Nhalf]) + np.min(lat[:Nhalf]))
            
        arr_c[npix, 0] = lon_c ; arr_c[npix, 1] = lat_c
        
    return arr_c

def find_points_inside_curve_(x, y, x_edge, y_edge, return_grid=False):
    """
    Find points enclosed by a curve.
    
    :param x = x-coordinate of points
    :param y = y-coordinate of points
    :param x_edge = x-coordinate of edge points (of curve/pixel)
    :param y_edge = y-coordinate of edge points (of curve/pixel)
    :param return_grid = boolean indicating output choice (default = False)
    
    :output 
        return_grid = False: (x,y) coordinates of points lying inside curve
        return_grid = True: indices (grid) of points that lie inside curve
        output.shape = (2, Ndata)
    
    Q: Why use return_grid = True?
    ANS: The "grid" is especially useful for thresholding the map of wavelet coefficients,
    allowing us to estimate wavelet coefficients only at designated regions. 
        (see generate_grid_points_ fct , generate_grid_points.py , and generate_wavelet_coefficients.py)
    """
    edge_points = np.stack((x_edge, y_edge), axis = -1) 
    edge_tuples = [tuple(edge_points[n,:]) for n in range(len(edge_points))]
    
    points = np.stack((x, y), axis = -1) 
    tuples = [tuple(points[n,:]) for n in range(len(points))]
    
    p = Path(edge_tuples) # make a polygon
    grid = p.contains_points(points)
    
    x_grid = x[grid]
    y_grid = y[grid]
    
    if return_grid == True:
        return grid
    else:
        return [x_grid, y_grid]

def find_points_inside_pixel_(x, y, npix, lon, lat, NSIDE=4, step=100):
    """
    Find data lying inside single father pixel. 
    
    :param x = longitude of data (rad)
    :param y = latitude of data (rad)
    :param npix = pixel number
    :param lon = longitude of father pixel edge
    :param lat = latitude of father pixel edge
    :param NSIDE = number of sides of healpix pixels (default = 4) 
    :param step = determines points needed to build pixel edges (default = 100)
    
    :output = list of longitude and latitude coordinates in father pixel
        [0][ndata] = list of longitudes indexed by npix (rad)
        [1][ndata] = list of latitudes indexed by npix (rad)
        output.shape = (2, Ndata)
    """

    disconnected_pixels = [40,72,104,136]
    lowest_pixels = [188,189,190,191]
    connected_right_pixels = [3,11,23,39,167,179,187,191]

    if npix not in disconnected_pixels:

        x_grid, y_grid = find_points_inside_curve_(x, y, lon, lat)

        ### Add points just below of lowest pixels to be counted
        if npix in lowest_pixels:
            n_x_miss = np.intersect1d(
            np.where(x > np.min(lon) )[0], np.where(x < np.max(lon))[0] )
            n_y_miss = np.where(y < np.min(lat))[0]
            n_miss = np.intersect1d(n_x_miss, n_y_miss)

            x_miss = x[n_miss]
            y_miss = y[n_miss]

            x_grid = np.concatenate((x_grid, x_miss))
            y_grid = np.concatenate((y_grid, y_miss))
        ### 

        ### Add points just to the right of rightmost, connected pixels to be counted
        if npix in connected_right_pixels:
            n_x_miss = np.where(x > np.max(lon))[0]
            n_y_miss = np.intersect1d(
            np.where(y > np.min(lat) )[0], np.where(y < np.max(lat))[0] )
            n_miss = np.intersect1d(n_x_miss, n_y_miss)

            x_miss = x[n_miss]
            y_miss = y[n_miss]

            x_grid = np.concatenate((x_grid, x_miss))
            y_grid = np.concatenate((y_grid, y_miss))
        ### 

    else:

        lon_noloop = np.delete(lon, -1)
        lat_noloop = np.delete(lat, -1)

        step = 100 # input for edge boundary creation
        N_edge_points = NSIDE * step # include extra point to close loop
        Nhalf = int(N_edge_points / 2)

        lon_1 = lon_noloop[0:Nhalf] # split data to two disjoint sets
        lon_2 = lon_noloop[Nhalf:]
        lat_1 = lat_noloop[0:Nhalf]
        lat_2 = lat_noloop[Nhalf:]

        lon_1 = np.append(lon_1, lon_2[0]) 
        lat_1 = np.append(lat_1, lat_2[0]) 
        lon_2 = np.delete(lon_2, 0) # remove line connecting _2 to _1
        lat_2 = np.delete(lat_2, 0)

        lon_1 = np.append(lon_1, lon_1[0]) # close up each region by completing loop
        lat_1 = np.append(lat_1, lat_1[0])
        lon_2 = np.append(lon_2, lon_2[0]) 
        lat_2 = np.append(lat_2, lat_2[0])

        x_grid_1, y_grid_1 = find_points_inside_curve_(x, y, lon_1, lat_1)
        x_grid_2, y_grid_2 = find_points_inside_curve_(x, y, lon_2, lat_2)

        x_grid = np.concatenate((x_grid_1, x_grid_2))
        y_grid = np.concatenate((y_grid_1, y_grid_2))

        ###
        # amendment: add points that are past the boundaries
        # of disjoint pixels
        n_x_miss_1 = np.where(x > np.max(lon_1))[0]
        n_y_between_1 = np.intersect1d(
            np.where(y > np.min(lat_1) )[0], np.where(y < np.max(lat_1))[0] )
        n_miss_1 = np.intersect1d(n_x_miss_1, n_y_between_1)
        x_miss_1 = x[n_miss_1]
        y_miss_1 = y[n_miss_1]

        n_x_miss_2 = np.where(x < np.min(lon_2))[0]
        n_y_between_2 = np.intersect1d(
            np.where(y > np.min(lat_2) )[0], np.where(y < np.max(lat_2))[0] )
        n_miss_2 = np.intersect1d(n_x_miss_2, n_y_between_2)
        x_miss_2 = x[n_miss_2]
        y_miss_2 = y[n_miss_2]

        ###

        x_grid = np.concatenate((x_grid, x_miss_1))
        x_grid = np.concatenate((x_grid, x_miss_2))
        y_grid = np.concatenate((y_grid, y_miss_1))
        y_grid = np.concatenate((y_grid, y_miss_2))
            
    return [x_grid, y_grid]

def divide_data_into_groups_(x, y, arr_edge_points, NSIDE=4, step=100):
    """
    Divide event data into groups based on which father pixel it occupies. 
    
    :param x = longitude of data (rad)
    :param y = latitude of data (rad)
    :param arr_edge_points = array of healpix edges (see healpix_edge_generator_ fct)
    :param NSIDE = number of sides of healpix pixels (default = 4) 
    :param step = determines points needed to build pixel edges (default = 100)
    
    :output = list of longitude and latitude coordinates labelled by father pixel
        [0][npix][ndata] = list of longitudes indexed by npix (rad)
        [1][npix][ndata] = list of latitudes indexed by npix (rad)
        output.shape = (2, Npix, Ndata)
    """
    NPIX, N_edge_points = arr_edge_points.shape[:-1] # number of father pixels, number of points in each pixel edge 
    
    grouped_points_x = []
    grouped_points_y = []

    disconnected_pixels = [40,72,104,136]
    lowest_pixels = [188,189,190,191]
    connected_right_pixels = [3,11,23,39,167,179,187,191]

    for npix in range(NPIX):
        if npix not in disconnected_pixels:
            lon = arr_edge_points[npix,:,0]
            lat = arr_edge_points[npix,:,1]

            x_grid, y_grid = find_points_inside_curve_(x, y, lon, lat)

            ### Add points just below of lowest pixels to be counted
            if npix in lowest_pixels:
                n_x_miss = np.intersect1d(
                np.where(x > np.min(lon) )[0], np.where(x < np.max(lon))[0] )
                n_y_miss = np.where(y < np.min(lat))[0]
                n_miss = np.intersect1d(n_x_miss, n_y_miss)

                x_miss = x[n_miss]
                y_miss = y[n_miss]

                x_grid = np.concatenate((x_grid, x_miss))
                y_grid = np.concatenate((y_grid, y_miss))
            ### 

            ### Add points just to the right of rightmost, connected pixels to be counted
            if npix in connected_right_pixels:
                n_x_miss = np.where(x > np.max(lon))[0]
                n_y_miss = np.intersect1d(
                np.where(y > np.min(lat) )[0], np.where(y < np.max(lat))[0] )
                n_miss = np.intersect1d(n_x_miss, n_y_miss)

                x_miss = x[n_miss]
                y_miss = y[n_miss]

                x_grid = np.concatenate((x_grid, x_miss))
                y_grid = np.concatenate((y_grid, y_miss))
            ### 

            grouped_points_x.append(x_grid)
            grouped_points_y.append(y_grid)

        else:
            lon = arr_edge_points[npix,:,0]
            lat = arr_edge_points[npix,:,1]

            lon_noloop = np.delete(lon, -1)
            lat_noloop = np.delete(lat, -1)

            step = 100 # input for edge boundary creation
            N_edge_points = NSIDE * step # include extra point to close loop
            Nhalf = int(N_edge_points / 2)

            lon_1 = lon_noloop[0:Nhalf] # split data to two disjoint sets
            lon_2 = lon_noloop[Nhalf:]
            lat_1 = lat_noloop[0:Nhalf]
            lat_2 = lat_noloop[Nhalf:]

            lon_1 = np.append(lon_1, lon_2[0]) 
            lat_1 = np.append(lat_1, lat_2[0]) 
            lon_2 = np.delete(lon_2, 0) # remove line connecting _2 to _1
            lat_2 = np.delete(lat_2, 0)

            lon_1 = np.append(lon_1, lon_1[0]) # close up each region by completing loop
            lat_1 = np.append(lat_1, lat_1[0])
            lon_2 = np.append(lon_2, lon_2[0]) 
            lat_2 = np.append(lat_2, lat_2[0])

            x_grid_1, y_grid_1 = find_points_inside_curve_(x, y, lon_1, lat_1)
            x_grid_2, y_grid_2 = find_points_inside_curve_(x, y, lon_2, lat_2)

            x_grid = np.concatenate((x_grid_1, x_grid_2))
            y_grid = np.concatenate((y_grid_1, y_grid_2))

            ###
            # amendment: add points that are past the boundaries
            # of disjoint pixels
            n_x_miss_1 = np.where(x > np.max(lon_1))[0]
            n_y_between_1 = np.intersect1d(
                np.where(y > np.min(lat_1) )[0], np.where(y < np.max(lat_1))[0] )
            n_miss_1 = np.intersect1d(n_x_miss_1, n_y_between_1)
            x_miss_1 = x[n_miss_1]
            y_miss_1 = y[n_miss_1]

            n_x_miss_2 = np.where(x < np.min(lon_2))[0]
            n_y_between_2 = np.intersect1d(
                np.where(y > np.min(lat_2) )[0], np.where(y < np.max(lat_2))[0] )
            n_miss_2 = np.intersect1d(n_x_miss_2, n_y_between_2)
            x_miss_2 = x[n_miss_2]
            y_miss_2 = y[n_miss_2]

            ###

            x_grid = np.concatenate((x_grid, x_miss_1))
            x_grid = np.concatenate((x_grid, x_miss_2))
            y_grid = np.concatenate((y_grid, y_miss_1))
            y_grid = np.concatenate((y_grid, y_miss_2))

            grouped_points_x.append(x_grid)
            grouped_points_y.append(y_grid)
            
    return [grouped_points_x, grouped_points_y]

def angular_separation(lon1, lat1, lon2, lat2): # taken from astropy.coordinates (later versions)
    """
    Angular separation between two points on a sphere.

    Parameters
    ----------
    lon1, lat1, lon2, lat2 : `~astropy.coordinates.Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the two points. Quantities should be in
        angular units; floats in radians.

    Returns
    -------
    angular separation : `~astropy.units.Quantity` ['angle'] or float
        Type depends on input; ``Quantity`` in angular units, or float in
        radians.

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1]_,
    which is slightly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    .. [1] https://en.wikipedia.org/wiki/Great-circle_distance
    """

    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.arctan2(np.hypot(num1, num2), denominator)

def find_neighboring_points_(ang_scale, x, y, x_edge, y_edge):
    """
    Finds point that are a set distance away from a curve (father pixel).
    
    :param ang_scale = angular distance threshold (rad)
    :param x = longitude of data (rad)
    :param y = latitude of data (rad)
    :param x_edge = longitude of edge points (of curve/pixel)
    :param y_edge = latitude of edge points (of curve/pixel)
    
    :output = (l,b)-coordinates of points that are ang_scale distance away from edges
        Note: The list will contain points that are also inside the curve
                  If undesirable, will need to remove these points with additional steps
                  (see remove_points_from_array_ as applied to generate_patches.py)
    """ 
    N_edge_points = len(x_edge)
    for n in range(N_edge_points):
        distances = angular_separation(x, y, x_edge[n], y_edge[n])
        x_out = x[distances < ang_scale]
        y_out = y[distances < ang_scale]
        if n == 0:
            points_out = np.stack((x_out, y_out), axis = -1)
        else:
            points_out = np.concatenate((points_out, np.stack((x_out, y_out), axis = -1)))

    # remove duplicates
    _, unique_indices = np.unique(points_out, axis = 0, return_index = True)
    x_out = points_out[unique_indices,0]
    y_out = points_out[unique_indices,1]
    
    return [x_out, y_out]

def remove_points_from_array_(A,B):
    """
    Ref: https://stackoverflow.com/questions/40055835/removing-elements-from-an-array-that-are-in-another-array
    
    Use the rows of B to remove rows of A (only tried for 2D arrays)
    """
    cumdims = (np.maximum(A.max(),B.max())+1)**np.arange(B.shape[1])
    return A[~np.in1d(A.dot(cumdims),B.dot(cumdims))]

def equal_ignore_order_(A, B):
    """
    Ref: https://stackoverflow.com/questions/40055835/removing-elements-from-an-array-that-are-in-another-array
    
    Show equality of arrays if they share same elements, independent of order (like set equality)
    """
    return np.array_equal(np.sort(A.flat), np.sort(B.flat))

def alpha_shape(points, alpha, only_outer=True):
    """
    * Ref: # https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points
    * alpha parameter reference: https://en.wikipedia.org/wiki/Alpha_shape
    
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n,2) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    
    """
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        """
        Add an edge between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                # if both neighboring triangles are in shape, it's not a boundary edge
                edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic = indices of corner points of the triangle
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        # Computing radius of triangle circumcircle
        # www.mathalino.com/reviewer/derivation-of-formulas/derivation-of-formula-for-radius-of-circumcircle
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return edges

def find_edges_with(i, edge_set):
    """
    * Ref: # https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points
    
    Necessary for stitch_boundaries file
    """
    i_first = [j for (x,j) in edge_set if x==i]
    i_second = [j for (j,x) in edge_set if x==i]
    return i_first,i_second

def stitch_boundaries(edges):
    """
    * Ref: # https://stackoverflow.com/questions/50549128/boundary-enclosing-a-given-set-of-points
    
    Stitches the output set of edge points from alpha_shape fct
    """
    edge_set = edges.copy()
    boundary_lst = []
    while len(edge_set) > 0:
        boundary = []
        edge0 = edge_set.pop()
        boundary.append(edge0)
        last_edge = edge0
        while len(edge_set) > 0:
            i,j = last_edge
            j_first, j_second = find_edges_with(j, edge_set)
            if j_first:
                edge_set.remove((j, j_first[0]))
                edge_with_j = (j, j_first[0])
                boundary.append(edge_with_j)
                last_edge = edge_with_j
            elif j_second:
                edge_set.remove((j_second[0], j))
                edge_with_j = (j, j_second[0])  # flip edge rep
                boundary.append(edge_with_j)
                last_edge = edge_with_j

            if edge0[0] == last_edge[1]:
                break

        boundary_lst.append(boundary)
    return boundary_lst

def generate_edge_of_point_set_(x, y, alpha=1):
    """
    Generates edge of a given point set using the "alpha-shape" formalism.
    
    :param x = x-coordinate of points
    :param y = y-coordinate of points
    :param alpha = parameter specifying the shape of the boundary (default = 1)
        (see alpha_shape fct)
        
    :output = array of points specifying the edge of the boundary set
    """
    points = np.stack((x, y), axis = -1)
    extension_edges = alpha_shape(points, alpha, only_outer = True)
    ext_edge = stitch_boundaries(extension_edges)[0]
    
    boundary_points = points[ext_edge,:][:,0,:]
    boundary_points = np.append(boundary_points, [boundary_points[0,:]], axis = 0)
    
    return boundary_points


# ======================================================================
# ======================================================================

# ======================================================================
# ======================================================================
                # MESH GRID GENERATION FUNCTIONS
    # References:
        # Development Directory: od/projects/pswavelets/gce/notebooks
            # grid_points_per_projected_map.ipynb
# ======================================================================
# ======================================================================

def build_mesh_(x_min, x_max, y_min, y_max, step_x, step_y, return_arrays_for_plotting=False):
    """
    Build 2D mesh given bounds and step sizes.
    
    :param x_min = x-grid lower bound 
    :param x_max = x-grid upper bound 
    :param y_min = y-grid lower bound 
    :param y_max = y-grid upper bound 
    :param step_x = x-grid step size
    :param step_y = y-grid step size
    :param return_arrays_for_plotting = boolean value for plotting
    
    :output 
        return_arrays_for_plotting = False: (x,y) mesh grid
        return_arrays_for_plotting = True: 
            [0] = (x,y) mesh grid
            [1] = array of (x,y) points making up grid
            [2] = array of lower-bound of "x-bins" used to define x-grid points
            [3] = array of lower-bound of "y-bins" used to define y-grid points
    """
    arr_x_plot = np.arange(x_min,x_max+step_x,step_x, dtype = float)
    arr_y_plot = np.arange(y_min,y_max+step_y,step_y, dtype = float)
    arr_x = 0.5 * (arr_x_plot[:-1] + arr_x_plot[1:])
    arr_y = 0.5 * (arr_y_plot[:-1] + arr_y_plot[1:])
    Nx = len(arr_x) ; Ny = len(arr_y)

    mesh_x, mesh_y = np.meshgrid(arr_x,arr_y) # each output array (NxN shaped) contains x or y value at given (i,j)-th position
    mesh_xy = np.stack((mesh_x, mesh_y), axis=-1)
    arr_r = mesh_xy.reshape(Nx*Ny,2) # flatten to 2D array
    
    if return_arrays_for_plotting == False:
        return mesh_xy
    else:
        return [mesh_xy, arr_r, arr_x_plot, arr_y_plot]
    
def generate_grid_points_(x_edge, y_edge, grid_scale, return_arrays_for_plotting=False):
    """
    Make grid points given the edge points of the outer region of a father pixel.
    
    :param x_edge = longitude of edge points (of curve/pixel)
    :param y_edge = latitude of edge points (of curve/pixel)
    :param grid_scale = step size of grid (step_x = step_y)
    :param return_arrays_for_plotting = boolean value for plotting
    
    :output
        return_arrays_for_plotting = False: 
            [0] = indices of grid points inside pixel edge (relative to mesh_bxby)
            [1] = (x,y) mesh grid
            [2] = list of points making up (x,y) mesh grid
       return_arrays_for_plotting = True:
           [3] = flattened set of indices given by [0]
           [4] = array of lower-bound of "x-bins" used to define x-grid points
           [5] = array of lower-bound of "y-bins" used to define y-grid points
    """
    bx_min, bx_max, by_min, by_max = [np.min(x_edge), np.max(x_edge),
                                      np.min(y_edge), np.max(y_edge)]
    step_size = grid_scale
    mesh_bxby, arr_b, arr_bx_plot, arr_by_plot = build_mesh_(bx_min, bx_max, by_min, by_max, step_size, step_size, True)
    Ny,Nx = mesh_bxby.shape[:-1]

#     # find grid points lying inside curve (useful only if node works only with a single grid point)
    grid_flat = find_points_inside_curve_(arr_b[:,0], arr_b[:,1], x_edge, y_edge, return_grid=True) # 2D array
    grid = grid_flat.reshape((Ny,Nx))
    if return_arrays_for_plotting == False:
        return [grid, mesh_bxby, arr_b]
    else:
        return [grid, mesh_bxby, arr_b, grid_flat, arr_bx_plot, arr_by_plot]

# ======================================================================
# ======================================================================
                # EXTRA FUNCTIONS
# ======================================================================
# ======================================================================

def closest_grid_point_(point, grid_points):
    relative_positions = grid_points - point
    distances = np.linalg.norm(relative_positions, axis = -1)
    bxby_index_shape = distances.shape
    flattened_index = np.argmin(distances)
    return np.unravel_index(flattened_index, bxby_index_shape)
    