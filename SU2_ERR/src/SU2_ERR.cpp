/*!
 * \file SU2_ERR.cpp
 * \brief Main file for the error computation code (SU2_ERR).
 * \author B. Munguía
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/SU2_ERR.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    
    unsigned short nDim, iZone, nZone = SINGLE_ZONE, iInst;
    unsigned long nPoint_coarse;
    su2double StartTime = 0.0, StopTime = 0.0, UsedTime = 0.0;
    ofstream ConvHist_file;
    char config_file_name[MAX_STRING_SIZE];
    int rank = MASTER_NODE;
    int size = SINGLE_NODE;
    bool periodic = false;
    
    /*--- MPI initialization, and buffer setting ---*/
    
#ifdef HAVE_MPI
    int  buffsize;
    char *buffptr;
    SU2_MPI::Init(&argc, &argv);
    SU2_MPI::Buffer_attach( malloc(BUFSIZE), BUFSIZE );
    SU2_Comm MPICommunicator(MPI_COMM_WORLD);
#else
    SU2_Comm MPICommunicator(0);
#endif
    rank = SU2_MPI::GetRank();
    size = SU2_MPI::GetSize();
    
    /*--- Pointer to different structures that will be used throughout the entire code ---*/
    
    CGeometry ****geometry_container        = NULL;
    CConfig **config_container              = NULL;
    CDriver *driver                         = NULL;
    unsigned short *nInst                  = NULL;
    
    /*--- Load in the number of zones and spatial dimensions in the mesh file (if no config
     file is specified, default.cfg is used) ---*/
    
    if (argc == 2 || argc == 3) { strcpy(config_file_name,argv[1]); }
    else { strcpy(config_file_name, "default.cfg"); }
    
    CConfig *config = NULL;
    config = new CConfig(config_file_name, SU2_ERR);
    
    nZone         = CConfig::GetnZone(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
    nDim          = CConfig::GetnDim(config->GetMesh_FileName(), config->GetMesh_FileFormat());
    periodic      = CConfig::GetPeriodic(config->GetMesh_FileName(), config->GetMesh_FileFormat(), config);
    nPoint_coarse = CConfig::GetnPoint(config->GetMesh_FileName(), config->GetMesh_FileFormat());
    
    /*--- Definition of the containers per zones ---*/
    
    config_container      = new CConfig*[nZone];
    geometry_container    = new CGeometry***[nZone];
    nInst = new unsigned short[nZone];
    
    for (iZone = 0; iZone < nZone; iZone++) {
        config_container[iZone]       = NULL;
        geometry_container[iZone]     = NULL;
        nInst[iZone]                  = 1;
    }
    
    /*--- Loop over all zones to initialize the various classes. In most
     cases, nZone is equal to one. This represents the solution of a partial
     differential equation on a single block, unstructured mesh. ---*/
    
    for (iZone = 0; iZone < nZone; iZone++) {
        
        /*--- Definition of the configuration option class for all zones. In this
         constructor, the input configuration file is parsed and all options are
         read and stored. ---*/
        
        config_container[iZone] = new CConfig(config_file_name, SU2_ERR, iZone, nZone, 0, VERB_HIGH);
        config_container[iZone]->SetMPICommunicator(MPICommunicator);
        
        /*--- Read the number of instances for each zone ---*/
        
        nInst[iZone] = config_container[iZone]->GetnTimeInstances();
        
        geometry_container[iZone] = new CGeometry**[nInst[iZone]];
        
        for (iInst = 0; iInst < nInst[iZone]; iInst++){
            
            /*--- Set some values for config. Note that for error estimation, we'll turn off MG. ---*/
            
            config_container[iZone]->SetiInst(iInst);
            config_container[iZone]->SetMGLevels(0);
            
            /*--- Definition of the geometry class to store the primal grid in the partitioning process. ---*/
            
            CGeometry *geometry_aux = NULL;
            
            /*--- All ranks process the grid and call ParMETIS for partitioning ---*/
            
            geometry_aux = new CPhysicalGeometry(config_container[iZone], iZone, nZone);
            
            /*--- Color the initial grid and set the send-receive domains (ParMETIS) ---*/
            
            geometry_aux->SetColorGrid_Parallel(config_container[iZone]);
            
            /*--- Allocate the memory of the current domain, and
             divide the grid between the nodes ---*/
            
            geometry_container[iZone][iInst] = NULL;
            
            geometry_container[iZone][iInst] = new CGeometry *[config_container[iZone]->GetnMGLevels()+1];
            
            /*--- Until we finish the new periodic BC implementation, use the old
             partitioning routines for cases with periodic BCs. The old routines
             will be entirely removed eventually in favor of the new methods. ---*/
            
            if (periodic) {
                geometry_container[iZone][iInst][MESH_0] = new CPhysicalGeometry(geometry_aux, config_container[iZone]);
            } else {
                geometry_container[iZone][iInst][MESH_0] = new CPhysicalGeometry(geometry_aux, config_container[iZone], periodic);
            }
            
            /*--- Deallocate the memory of geometry_aux ---*/
            
            delete geometry_aux;
            
            /*--- Add the Send/Receive boundaries ---*/
            
            geometry_container[iZone][iInst][MESH_0]->SetSendReceive(config_container[iZone]);
            
            /*--- Add the Send/Receive boundaries ---*/
            
            geometry_container[iZone][iInst][MESH_0]->SetBoundaries(config_container[iZone]);
            
            /*--- Create the vertex structure (required for MPI) ---*/
            
            if (rank == MASTER_NODE) cout << "Identify vertices." <<endl;
            geometry_container[iZone][iInst][MESH_0]->SetVertex(config_container[iZone]);
            
            /*--- Store the global to local mapping after preprocessing. ---*/
            
            if (rank == MASTER_NODE) cout << "Storing a mapping from global to local point index." << endl;
            geometry_container[iZone][iInst][MESH_0]->SetGlobal_to_Local_Point();
            
        }
        
    }
    
    delete config;
    config = NULL;
    
    unsigned long nPoint_fine = geometry_container[ZONE_0][INST_0][MESH_0]->GetGlobal_nPointDomain();
    if(rank == MASTER_NODE) cout << "Points in coarse mesh = " << nPoint_coarse << "\nPoints in fine mesh   = " << nPoint_fine << endl;
    
    /*--- Set up a timer for performance benchmarking (preprocessing time is included) ---*/
    
#ifdef HAVE_MPI
    StartTime = MPI_Wtime();
#else
    StartTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif
    
    /*--- Define a discrete adjoint driver with the refined geometry and updated config options. ---*/
    
    driver = new CDiscAdjFluidDriver(config_container, geometry_container, nZone, nDim, periodic, MPICommunicator);
    
    /*--- Perform the error estimation. ---*/
    driver->ComputeErrorEstimate(nPoint_coarse, nPoint_fine);
    
    delete driver;
    
    if (rank == MASTER_NODE)
        cout << endl <<"------------------------- Solver Postprocessing -------------------------" << endl;
    
    if (geometry_container != NULL) {
        for (iZone = 0; iZone < nZone; iZone++) {
            for (iInst = 0; iInst < nInst[iZone]; iInst++){
                if (geometry_container[iZone][iInst] != NULL) {
                    delete geometry_container[iZone][iInst];
                }
            }
            if (geometry_container[iZone] != NULL)
                delete geometry_container[iZone];
        }
        delete [] geometry_container;
    }
    if (rank == MASTER_NODE) cout << "Deleted CGeometry container." << endl;
    
    if (config_container != NULL) {
        for (iZone = 0; iZone < nZone; iZone++) {
            if (config_container[iZone] != NULL) {
                delete config_container[iZone];
            }
        }
        delete [] config_container;
    }
    if (rank == MASTER_NODE) cout << "Deleted CConfig container." << endl;
    
    /*--- Synchronization point after a single solver iteration. Compute the
     wall clock time required. ---*/
    
#ifdef HAVE_MPI
    StopTime = MPI_Wtime();
#else
    StopTime = su2double(clock())/su2double(CLOCKS_PER_SEC);
#endif
    
    /*--- Compute/print the total time for performance benchmarking. ---*/
    
    UsedTime = StopTime-StartTime;
    if (rank == MASTER_NODE) {
        cout << "\nCompleted in " << fixed << UsedTime << " seconds on "<< size;
        if (size == 1) cout << " core." << endl; else cout << " cores." << endl;
    }
    
    /*--- Exit the solver cleanly ---*/
    
    if (rank == MASTER_NODE)
        cout << endl <<"------------------------- Exit Success (SU2_ERR) ------------------------" << endl << endl;
    
    /*--- Finalize MPI parallelization ---*/
    
#ifdef HAVE_MPI
    SU2_MPI::Finalize();
#endif
    
    return EXIT_SUCCESS;
}
