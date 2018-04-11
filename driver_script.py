''' This script generates  5 RVEs of plain woven composites with 11%, 18%, 32%,
41% and 53% yarn volume fraction.
'''
import os
import sys
import shutil
import numpy as np
import mesh_lib
################################################################################
# INPUT PARAMETERS
################################################################################
# ABAQUS-PATH
abq_path ='/opt/SIMULIA/Commands/abaqus '
# INPUT-FILE-NAME
inp_file_name_tmp = 'RVE_plain_woven_PMESH'
################################################################################
# ADUJSTABLE GEOMETRY PARAMETERS
v_fracs        =              [11   ,18  ,32   ,41  ,53  ]
ellipse_widths =              [0.2  ,0.32,0.55 ,0.7 ,0.9 ] # width of ellipse
# MESHING PARAMETERS
mesh_densities_matrix       = [0.09 ,0.07,0.07 ,0.09,0.09 ]
mesh_densities_yarn_ellipse = [0.055,0.06,0.035,0.05,0.05 ]
mesh_densities_yarn_inside  = [0.065,0.08,0.04 ,0.06,0.06 ]
quad_ele    = True
ele_ord_str = 'QUAD'
# MATERIAL DATA
isotropy_flag   = True
isotropy_string = 'ISO'
E_Matrix = 4350.
nu_Matrix= 0.36
E_Yarn   = 10590.
nu_Yarn  = 0.15
# LOADING-DATA
epsilon = 0.02 # strain for macroscopic elasticity tensor
u       = 0.1  # displacement of tensile test

################################################################################
################################################################################
# PROGRAM EXECUTION
################################################################################
################################################################################
cmd     = abq_path+' cae noGUI=gen_plain_woven_RVE_in_abq.py -- "'
print 'ABAQUS-PATH: '+abq_path
print '========================================================================'
# GENERATE INPUT FILES
for v_frac,b,mesh_density_matrix,mesh_density_yarn_ellipse,mesh_density_yarn_inside in zip(v_fracs,ellipse_widths,mesh_densities_matrix,mesh_densities_yarn_ellipse,mesh_densities_yarn_inside):
    print 'processing inclusion volume fraction:', v_frac, '%'
    inp_file_name = inp_file_name_tmp+'_vfrac'+str(v_frac)
    # check if inp-file exists allready
    continue_flag = False
    for root, directories, filenames in os.walk(os.getcwd()):
        for filename in filenames:
            if filename.endswith('.inp'):
                if filename.startswith(inp_file_name):
                    continue_flag=True
    if continue_flag:
        continue
    # yarn cross-section        
    yrn_cmd  = ' b='+str(b)+' ' # width of ellipse    
    yrn_cmd += ' h=0.4'+' '     # height of ellipse
    yrn_cmd += ' d1=0.06'+' '   # yarn distance in y-direction at maximum
    yrn_cmd += ' T=2.5'+' '     # period length
    # meshing_parameters
    msh_cmd  = ' mesh_density_matrix='      +str(mesh_density_matrix)      +' '
    msh_cmd += ' mesh_density_yarn_ellipse='+str(mesh_density_yarn_ellipse)+' '
    msh_cmd += ' mesh_density_yarn_inside=' +str(mesh_density_yarn_inside) +' '
    # GENERATE BASIC INPUT FILE
    final_cmd = cmd+msh_cmd+' inp_file_name='+inp_file_name+' '+yrn_cmd+' "'    
    if inp_file_name+'.inp' not in os.listdir(os.getcwd()):
        os.system(final_cmd)
        mesh = mesh_lib.read_abq_inp(inp_file_name)    
        part = mesh.parts[mesh.parts.keys()[0]]
        # take care of abq-counting
        for key in part.ele_sets.keys():
            part.ele_sets[key] = list(np.array(part.ele_sets[key])-1)
        mesh_lib.mesh2abq(inp_file_name, part.nds, part.elements['C3D4']['conn'],ele_sets=part.ele_sets)
    #CONVERT MESH TO QUADRATIC ELEMENTS
    if quad_ele:
        print '    convert to quadratic elements'
        ele_ord_str = 'QUAD'
        inp_file_name_final = inp_file_name+'_'+ele_ord_str+'TET'
        mesh_lib.conv_lin_abq_inp_msh2_quad_abq_inp(inp_file_name,inp_file_name_final)
    else:
        ele_ord_str = 'LIN'
        inp_file_name_final = inp_file_name+'_'+ele_ord_str++'TET'
        shutil.copy(inp_file_name+'.inp', inp_file_name_final+'.inp')
    ################################################################################
    print 'applying periodic boundary conditions'
    # APPLY PERIODIC BOUNDARY CONDITIONS    
    mesh_lib.apply_PBC(inp_file_name_final,u,epsilon,
                       isotropy_flag=True,
                       E_Matrix=4350.,nu_Matrix=0.36,
                       E_Yarn  =10590.,nu_Yarn=0.15)
    print '----------------'