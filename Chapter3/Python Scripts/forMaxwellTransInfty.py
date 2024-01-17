
import time
import os.path
t = time.time()

import numpy as np
import bempp.api
import scipy.io as sio
bempp.api.enable_console_logging()
# BEMPP Parameters
bempp.api.global_parameters.assembly.potential_operator_assembly_type = 'dense'

radii = ['003','004','005','006','007','008','009','01','012','014','016','018','02']
#
filename = 'runtime_information_infinity.txt'
with open(filename, 'w') as f:
    f.write('Script about the runtime information\n')
bempp.api.global_parameters.hmat.eps = 1E-9
frequency = 100E6
EpsExt = 8.854187817E-12
MuExt = 4 * np.pi * 1E-7
EpsRel = 1.0
MuRel = 2.1
WavenumberExt = 2 * np.pi * frequency * np.sqrt(MuExt * EpsExt)
WavenumberInt = WavenumberExt * np.sqrt(EpsRel*MuRel)
Points = np.loadtxt('example.txt')
direction = np.array([1/np.sqrt(3),-1/np.sqrt(3),1/np.sqrt(3)])
polarization = np.array([-1.0, 1.0j, 1.0+1.0j])
with open(filename, 'a') as f:
    f.write('hmateps = ' + str(1E-9) + '\n' +
            'Eps = ' + str(EpsExt) + '\n' +
            'Mu = ' + str(MuExt) + '\n' +
            'Epsrel = ' + str(EpsRel) + '\n' +
            'Murel = ' + str(MuRel) + '\n' +
            'WavenumberExt = ' + str(WavenumberExt) + '\n' +
            'WavenumberInt = ' + str(WavenumberInt) + '\n' +  
            'direction = ' + str(direction) + '\n' +
            'polarization = ' + str(polarization) + '\n' 
            )
# 
# GMRES information
tol = 1e-5
restart = None
use_strong_form = True
return_iteration_count = True
for x in radii:
    name_of_ff = os.path.join('farfields/dielectric/infinity', 'farFieldinfty_new'+x+'.mat')
    name_of_mesh = os.path.join('grids/infinity', 'infty_new' + x + '.msh')
    ##
    with open(filename, 'a') as f:
        f.write('Radius r = ' + str(x) + '\n')
    
    def plane_wave(point):
        return polarization * np.exp(1j * WavenumberExt * np.dot(point, direction))
	
    def scaled_plane_wave(point):
        return np.sqrt(EpsExt) * plane_wave(point)
	
    def tangential_trace(point, n, domain_index, result):
        result[:] =  np.cross(scaled_plane_wave(point), n)
	
    def scaled_plane_wave_curl(point):
        return np.cross(direction, polarization) * 1j * WavenumberExt * np.sqrt(EpsExt) * np.exp(1j * WavenumberExt * np.dot(point, direction))
	
    def neumann_trace(point, n, domain_index, result):
        result[:] =  1./ (1j * WavenumberExt) * np.cross(scaled_plane_wave_curl(point), n)
	
    start = time.time()
    Grid=bempp.api.import_grid(name_of_mesh)
    end = time.time()
    elapsed_time = end-start
    with open(filename, 'a') as f:
        f.write('It took ' + str(elapsed_time) + ' seconds to load the grid \n')

    start = time.time()
    MultitraceInt = bempp.api.operators.boundary.maxwell.multitrace_operator(Grid, WavenumberInt)
    MultitraceExt = bempp.api.operators.boundary.maxwell.multitrace_operator(Grid, WavenumberExt)
	
    def ScalingSASinv(A):
        A[0,1] = A[0,1] * np.sqrt(MuRel / EpsRel)
        A[1,0] = A[1,0] * np.sqrt(EpsRel / MuRel)
        return A
	
    ScaledMultitraceInt = ScalingSASinv(MultitraceInt)
    OperatorLHS = ScaledMultitraceInt + MultitraceExt
    end = time.time()
    elapsed_time = end-start
    with open(filename, 'a') as f:
        f.write('It took ' + str(elapsed_time) + ' seconds to set up the left hand side of the system \n')

    start = time.time()
    Identity = bempp.api.operators.boundary.sparse.multitrace_identity(Grid, spaces='maxwell')
    OperatorRHS = 0.5 * Identity - ScaledMultitraceInt
    end = time.time()
    elapsed_time = end-start
    with open(filename, 'a') as f:
        f.write('It took ' + str(elapsed_time) + ' seconds to set up the right hand side of the system with scaledmultitraceint given \n')
    IncidentTraces = (bempp.api.GridFunction(OperatorRHS.domain_spaces[0], fun=tangential_trace, dual_space = OperatorRHS.dual_to_range_spaces[0]),
	                  bempp.api.GridFunction(OperatorRHS.domain_spaces[1], fun=neumann_trace, dual_space = OperatorRHS.dual_to_range_spaces[1]))
    Rhs = OperatorRHS * IncidentTraces
    start = time.time()
    Solution, info, numofit = bempp.api.linalg.gmres(OperatorLHS * OperatorLHS, OperatorLHS * Rhs, tol=tol, restart =restart,
                                            return_iteration_count=return_iteration_count, use_strong_form=use_strong_form)

    end = time.time()
    elapsed_time = end-start
    with open(filename, 'a') as f:
        f.write('Parameters of GMRES:\n'+
                'Tolerance = '+ str(tol)+'\n'+
                'Restart = '+ str(restart)+'\n'+
                'Strong form = '+ str(use_strong_form)+'\n'+
                'Number of dofs = ' + str(OperatorLHS.domain_spaces[0].global_dof_count)+'\n'+
                'Number of iterations = '+ str(numofit)+'\n'+
                'It took ' + str(round(elapsed_time,2)) + ' seconds to perform GMRES\n')
    FarFieldOperator = (bempp.api.operators.far_field.maxwell.magnetic_field(Solution[0].space, Points, WavenumberExt),
                        bempp.api.operators.far_field.maxwell.electric_field(Solution[1].space, Points, WavenumberExt))
    FarField = - 4.0 * np.pi * (FarFieldOperator[0] * Solution[0] + FarFieldOperator[1] * Solution[1])
    sio.savemat(name_of_ff,{'FarField':FarField})