# ===================================================================
# This is a user configuration file for ParMooN Version 1.1
# written by Sashikumaar Ganesan, CDS, IISc Bangalore, India
# date: 05 June 2015
# ===================================================================

# controlling the output messages
set(CMAKE_VERBOSE_MAKEFILE FALSE)

# selection of dimension (2D 3D)
set(AParMooN_GEO "2D" CACHE STRING "Change AParMooN_GEO, to select the Dimensio of the problem")
# set(AParMooN_GEO "2D" CACHE STRING "Change AParMooN_GEO, to select the Dimension of the problem")

#...................................................................................................................................................
####################################------Selection of main program start ------##################################
# Uncomment the main program you would like to run:

##### -------------------------------------------(A) DO Programs------------------------------------------------#####
#1. Time Dependent Linear Advection DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/TLA2D_ParMooN_IVUQ_DO.C" CACHE STRING "Enter to select the Main file of the model") 

#2. Time Dependent Convection-Diffusion DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/TCD2D_ParMooN_IVUQ_DO.C" CACHE STRING "Enter to select the Main file of the model")

#3. Time Dependent Burgers' DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/TBE2D_ParMooN_IVUQ_DO.C" CACHE STRING "Enter to select the Main file of the model")

#4. Time Dependent Navier-Stokes' DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/TNSE2D_ParMooN_IVUQ_DO.C" CACHE STRING "Enter to select the Main file of the model")


##### -----------------------------------------(B) DO Test Programs----------------------------------------------#####
#1. Time Dependent Linear Advection DO 
set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/Trials/TLA2D_ParMooN_IVUQ_DO_Test.C" CACHE STRING "Enter to select the Main file of the model") 

#2. Time Dependent Convection-Diffusion DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/Trials/TCD2D_ParMooN_IVUQ_DO_Test.C" CACHE STRING "Enter to select the Main file of the model")

#3. Time Dependent Burgers' DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/Trials/TBE2D_ParMooN_IVUQ_DO_Test.C" CACHE STRING "Enter to select the Main file of the model")

#4. Time Dependent Navier-Stokes' DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/Trials/TNSE2D_ParMooN_IVUQ_DO_Test.C" CACHE STRING "Enter to select the Main file of the model")



#5. Miscellaneous test programs
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/TCD2D_ParMooN_IVUQ_DO_Mean.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/TBE2D_ParMooN_IVUQ_DO_Mean.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TBE2D_ParMooN_MC_Test.C" CACHE STRING "Enter to select the Main file of the model") 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/DO_UQ/TNSE2D_ParMooN_IVUQ_DO_Test.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/TESTS/TNSE2D_ParMooN_IVUQ_DO_NoDF.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/TESTS/TNSE2D_ParMooN_IVUQ_DO_DFPre.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/TESTS/TNSE2D_ParMooN_IVUQ_DO_DFPost.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/TESTS/TNSE2D_ParMooN_IVUQ_MC_wMean.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/TESTS/TCD2D_ParMooN_MC_Mean.C" CACHE STRING "Enter to select the Main file of the model")

##### ---------------------------------------(C) Monte Carlo Programs--------------------------------------------#####
#1. Time Dependent Linear Advection DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TLA2D_ParMooN_MC.C" CACHE STRING "Enter to select the Main file of the model")

#2. Time Dependent Convection-Diffusion DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TCD2D_ParMooN_MC.C" CACHE STRING "Enter to select the Main file of the model")

#3. Time Dependent Burgers' DO 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TBE2D_ParMooN_MC.C" CACHE STRING "Enter to select the Main file of the model") 

#4. Time Dependent Navier-Stokes' DO (Initial Value Uncertainty)
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TNSE2D_ParMooN_IVUQ_MC.C" CACHE STRING "Enter to select the Main file of the model")
#4. Time Dependent Navier-Stokes' DO (Parameter Uncertainty)
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TNSE2D_ParMooN_RE_NR_MC.C" CACHE STRING "Enter to select the Main file of the model")


### 2 Dimensional Problems ###

#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/CD2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model") #Steady State Convection-Diffusion Equation

#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/TNSE2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model") #Time-dependent Navier-Stokes' Equation

# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/TCD2D_ParMooN_ALE.C" CACHE STRING "Enter to select the Main file of the model") 
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/NSE2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")#Steady State Navier-Stokes' Equation
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/TNSE2D_IVUQ_SPADE.C" CACHE STRING "Enter to select the Main file of the model")

#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/TBE2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model") #Time-dependent Burgers' Equation
### 3 Dimensional Problems ###
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/CD3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/TCD3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/NSE3D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/3DPrograms/TNSE2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model")


# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Main_Users/Thivin/2D_Programs/deepikaji.cpp" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Main_Users/Thivin/TNSE3D/thivin_TNSE3D.cpp" CACHE STRING "Enter to select the Main file of the model")
# set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Main_Users/Thivin/Sample_mesh_move_2d.cpp" CACHE STRING "Enter to select the Main file of the model")

#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/2DPrograms/TBE2D_ParMooN.C" CACHE STRING "Enter to select the Main file of the model") 
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/1DPrograms/CD1D_ParMooN_ANN.C" CACHE STRING "Enter to select the Main file of the model") 

##### Monte Carlo Programs#####

#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TCD2D_ParMooN_MC.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TNSE2D_ParMooN_IVUQ_MC.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TNSE2D_ParMooN_RE_NR_MC.C" CACHE STRING "Enter to select the Main file of the model")
#set(AParMooN_MODEL "${PROJECT_SOURCE_DIR}/Monte_Carlo/TBE2D_ParMooN_MC.C" CACHE STRING "Enter to select the Main file of the model") #Time-dependent Burgers' Equation






############################### ---Selection of main program end ----############################################



# selection of architect type (LINUX64 MAC64 INTEL64 TYRONE64 CRAY64)
set(AParMooN_ARCH "LINUX64" CACHE STRING "select the machine type")

#  selection of program type (SEQUENTIAL SMPI MPI OMPONLY HYBRID SCUDA)
set(AParMooN_PARALLEL_TYPE "SEQUENTIAL" CACHE STRING "select the parallel type")

#  selection of program type (MPICH OPENMPI INTELMPI CRAYMPI MACMPI)
set(AParMooN_MPI_IMPLEMENTATION "INTELMPI" CACHE STRING "select the MPI Implementation type")


############################### Selection of output folder start ################################################
# set the path to save the exe file ....................................................................................
#.......................................................................................................................
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/aletnse3d" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/TCD2D_MC" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/NSE2D" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/TNSE2D" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/CD2D" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/TCD3D" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/NSE3D" CACHE STRING "select the model")

# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/TNSE2D" CACHE STRING "select the model")

# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/thivin3d" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/mesh2d" CACHE STRING "select the model")
# set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/deepikaji" CACHE STRING "select the model")

#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/TBE2D" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/cd1dANN" CACHE STRING "select the model")

####### Monte-Carlo Output Folders ##########
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Monte_Carlo/TCD2D_MC" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Monte_Carlo/TBE2D_MC" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Monte_Carlo/TNSE2D_IVUQ_MC" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Monte_Carlo/TNSE2D_PUQ_MC" CACHE STRING "select the model")

####### DO Output Folders ##########
#1. Output folder for time-dependent linear advection DO solution
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TLA2D_IVUQ_DO" CACHE STRING "select the model") 

#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TNSE2D_IVUQ_DO" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TCD2D_IVUQ_DO" CACHE STRING "select the model")


##### Test Output Folders ###########################
set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TLA2D_IVUQ_DO_Test_New" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TCD2D_IVUQ_DO_Test1" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TCD2D_IVUQ_DO_Mean" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TBE2D_IVUQ_DO_Test" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TBE2D_IVUQ_DO_Mean" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Monte_Carlo/TBE2D_MC_Test" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/DO/TNSE2D_IVUQ_DO_Test" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Test/TNSE2D_IVUQ_DO_NoDF" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Test/TNSE2D_IVUQ_DO_DFPre" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Test/TNSE2D_IVUQ_DO_DFPost" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Test/TNSE2D_IVUQ_MC" CACHE STRING "select the model")
#set(AParMooN_OUTPUT_DIR_PATH "${CMAKE_SOURCE_DIR}/../ParMooN_Output/Test/TCD2D_MC_Mean" CACHE STRING "select the model")



############################### Selection of output folder end ################################################
#--------------------------------------------------------------------------------------------------------------------
set(USE_PARMOON_DEFINE -D__PRIVATE__)

#------------------------------------------------------------------------------------------------------------------------
# CMAKE_BUILD_TYPE [ DEBUG | RELEASE | RELWITHDEBINFO | MINSIZEREL ]
set(EXE_BUILD_TYPE DEBUG)


# set FALSE, if you want to use libs provided in PARMOON 
# if you set TRUE, it will search in all you lib paths and if not found, PARMOON libs will be used
set(AParMooN_USE_SYSTEM_MKLBLAS TRUE)
set(AParMooN_USE_SYSTEM_UMFPACK TRUE)
set(AParMooN_USE_SYSTEM_LAPACK TRUE)
set(AParMooN_USE_SYSTEM_MUMPS TRUE)
set(AParMooN_USE_SYSTEM_GRIDGEN TRUE)
set(AParMooN_USE_SYSTEM_TETGEN TRUE)

# set(USE_PARMOON_DEFINE -D__PRIVATE__)
# SET(DAI_WITH_BP ON CACHE BOOL "Belief Propagation" FORCE)
# ========================================================================================================================
# no need to change anyting after this line
# used only when ccmake or cmake-gui is used
# ========================================================================================================================
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${AParMooN_OUTPUT_DIR_PATH}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${AParMooN_OUTPUT_DIR_PATH}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${AParMooN_OUTPUT_DIR_PATH})

set_property(CACHE AParMooN_GEO PROPERTY STRINGS 2D 3D ) 

# selection of all main programs
if("${AParMooN_GEO}" STREQUAL "2D")
  file(GLOB_RECURSE MAIN_SOURCES "${PROJECT_SOURCE_DIR}/2DPrograms/*.C")
  set_property(CACHE AParMooN_MODEL PROPERTY STRINGS  ${MAIN_SOURCES})   
elseif("${AParMooN_GEO}" STREQUAL "3D")
  file(GLOB_RECURSE MAIN_SOURCES "${PROJECT_SOURCE_DIR}/3DPrograms/*.C")
  set_property(CACHE AParMooN_MODEL PROPERTY STRINGS  ${MAIN_SOURCES})  
endif()  

# selection of all architects
set_property(CACHE AParMooN_ARCH PROPERTY STRINGS LINUX64 MAC64 INTEL64 TYRONE64 CRAY64)

# selection of all program types
set_property(CACHE AParMooN_PARALLEL_TYPE PROPERTY STRINGS SEQUENTIAL SMPI MPI OMPONLY HYBRID SCUDA)

# selection of all program types
set_property(CACHE AParMooN_MPI_IMPLEMENTATION PROPERTY STRINGS MPICH OPENMPI INTELMPI CRAYMPI MACMPI)

# selection of all output folder
if("${AParMooN_GEO}" STREQUAL "2D")
set_property(CACHE AParMooN_OUTPUT_DIR_PATH PROPERTY STRINGS  ${CMAKE_SOURCE_DIR}/OutPut/cd2d  ${CMAKE_SOURCE_DIR}/OutPut/tcd2d  ${CMAKE_SOURCE_DIR}/OutPut/nse2d  ${CMAKE_SOURCE_DIR}/OutPut/tnse2d)
 elseif("${AParMooN_GEO}" STREQUAL "3D")
set_property(CACHE AParMooN_OUTPUT_DIR_PATH PROPERTY STRINGS  ${CMAKE_SOURCE_DIR}/OutPut/cd3d  ${CMAKE_SOURCE_DIR}/OutPut/tcd3d  ${CMAKE_SOURCE_DIR}/OutPut/nse3d  ${CMAKE_SOURCE_DIR}/OutPut/tnse3d)
endif()

# ======================================================================
# general settings
# ======================================================================
 if("${AParMooN_PARALLEL_TYPE}" STREQUAL "MPI")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_MPIONLY -D_MPI")
 elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "OMPONLY")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_OMPONLY -D_OMP")
 elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "HYBRID")
   set(PARMOON_PRG_DEFINE "-D_PAR -D_MPI -D_HYBRID")
 elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SEQUENTIAL")
    set(PARMOON_PRG_DEFINE "-D_SEQ ")
 elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SCUDA")
    set(PARMOON_PRG_DEFINE "-D_SEQ -D_CUDA")
  elseif("${AParMooN_PARALLEL_TYPE}" STREQUAL "SMPI")
    set(PARMOON_PRG_DEFINE "-D_SEQ -D_SMPI")
 endif()

 if("${AParMooN_ARCH}" STREQUAL "LINUX64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -fopenmp -std=c++11 -lmkl_rt -lmkl_intel_thread -lpthread -liomp5")
   set(PARMOON_C_DEF "  ${PARMOONx_C_DEF}  -DREDUCED -DNO_TIMER -m64  -lmkl_rt -lmkl_intel_thread  ")
   set(PARMOON_CUDA_DEF "${PARMOON_CUDA_DEF} -arch=sm_70 -Xcompiler -fopenmp -w ")
   
 elseif("${AParMooN_ARCH}" STREQUAL "MAC64")
   set(PARMOON_CXX_DEF " -stdlib=libc++ ${PARMOON_CXX_DEF}  -fapple-pragma-pack -Wdeprecated-register  ")
   set(PARMOON_C_DEF "  ${PARMOON_C_DEF}  -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64 -fapple-pragma-pack ")
 elseif("${AParMooN_ARCH}" STREQUAL "INTEL64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -std=c++11   -qopenmp ")
   set(PARMOON_C_DEF "  ${PARMOON_C_DEF}   -DREDUCED -DNO_TIMER -DMKL_ILP64 -m64  ")
   set(PARMOON_CUDA_DEF "${PARMOON_CUDA_DEF} -arch=sm_70 -Xcompiler -fopenmp -w ")

 elseif("${AParMooN_ARCH}" STREQUAL "TYRONE64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -DREDUCED -DNO_TIMER")
   set(PARMOON_C_DEF "  ${PARMOON_C_DEF}  -DREDUCED -DNO_TIMER  -DMPICH_IGNORE_CXX_SEEK ")  
 elseif("${AParMooN_ARCH}" STREQUAL "CRAY64")
   set(PARMOON_CXX_DEF "${PARMOON_CXX_DEF} -DMPICH_IGNORE_CXX_SEEK ")  
   set(PARMOON_C_DEF "  ${PARMOON_C_DEF}  -DREDUCED -DNO_TIMER -DMPICH_IGNORE_CXX_SEEK  ")
   set(PARMOON_CUDA_DEF "${PARMOON_CUDA_DEF} -arch=sm_70 -Xcompiler -fopenmp -w ")

 endif()
 
set(PARMOON_C_DEF " ${PARMOON_C_DEF} -D__${AParMooN_GEO}__ -D__${AParMooN_ARCH}__  -DTRILIBRARY -DTETLIBRARY  ${PARMOON_PRG_DEFINE} ${USE_PARMOON_DEFINE}") 
set(PARMOON_CXX_DEF " ${PARMOON_CXX_DEF} -D__${AParMooN_GEO}__ -D__${AParMooN_ARCH}__  -DTRILIBRARY -DTETLIBRARY ${PARMOON_PRG_DEFINE} ${USE_PARMOON_DEFINE}")
set(PARMOON_CUDA_DEF " ${PARMOON_CUDA_DEF} -D__${AParMooN_GEO}__ -D__${AParMooN_ARCH}__  -DTRILIBRARY -DTETLIBRARY ${PARMOON_PRG_DEFINE} ${USE_PARMOON_DEFINE}")
