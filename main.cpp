/* *********************************************************************
 * Universidade Federal Fluminense (UFF)
 * 
 * The code is build usind PETSc, for more information acess:
 * https://www.mcs.anl.gov/petsc/
 * 
 * One can use and modify this code for educational purpose only.
 * 
 * Created By: Ricardo Dias dos Santos <ricardos@id.uff.br>
 * Date: 09/26/2018
 **********************************************************************/

#include <petscsys.h> 
#include <petsctime.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmlabel.h> 
#include <petscds.h>
#include <petscdraw.h>
#include <petscts.h>

#define _MAX_STRING_SIZE_ 200

template <typename type>
PETSC_STATIC_INLINE type pow2(type valor){
    return valor*valor;
}

template <typename type>
PETSC_STATIC_INLINE type pow3(type valor){
    return valor*valor*valor;
}

static TSARKIMEXType  TSSDIRK2    = "SDIRK22";
static TSARKIMEXType  TSSDIRK3    = "SDIRK33";
static TSARKIMEXType  TSSDIRK2SSP = "SDIRK22SSP";
static TSARKIMEXType  TSSDIRK3SSP = "SDIRK33SSP";
static TSARKIMEXType  TSRKMGM2    = "SDIRK2MGM";
static TSARKIMEXType  TSRKMGM3    = "SDIRK3MGM";
static PetscErrorCode RegisterMyRK(PetscReal);

//#include <stdio.h>
//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <cstring>
#include <string>

#include "lib/data_strucuture/data_strucuture.hpp"
#include "lib/parameters/parameters.hpp"
#include "lib/dmda/dmda.hpp"
#include "lib/vectors/vectors.hpp"
#include "lib/simulation_parameters/simulation_parameters.hpp" 
#include "lib/linear_algebra/linear_algebra.hpp"
#include "lib/boostconv/boostconv.hpp"
#include "lib/interpolation/interpolation.hpp"
#include "lib/initial_condition/initial_condition.hpp"
#include "lib/boundary_condition/boundary_condition.hpp"
#include "lib/entropy_fix/entropy_fix.hpp"
#include "lib/reconstruction/reconstruction.hpp"
#include "lib/flux/flux.hpp"
#include "lib/flux_splitting/flux_splitting.hpp"
#include "lib/jacobian/jacobian.hpp"
#include "lib/derivative/derivative.hpp"
#include "lib/viscous/viscous.hpp"
#include "lib/buffer/buffer.hpp"
#include "lib/perturbation/perturbation.hpp"
#include "lib/source/source.hpp"
#include "lib/rhs/rhs.hpp"
#include "lib/output/output.hpp"
#include "lib/time_integration/time_integration.hpp"
#include "lib/nan/nan.hpp"
#include "lib/initialize/initialize.hpp"

static char help[] = "Solves a Compressible 1D flow.\n\n";

int main(int argc,char **args){
    
    PetscErrorCode  ierr;
    PetscLogDouble  t1,t2;
    Solve           *solve;
    Functions       *func;
    DAVec           *dmda;
    Buffer			*buffer;
    Pert            *pert;
    Viewers			*view;
    Flags			*flags;
    
    PetscInitialize(&argc,&args,(char*)0,help);
    
    /*******************************************************************
     *                   Allocate Solver Structs
     ******************************************************************/
    
    ierr = PetscMalloc1(1,&solve);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->prmt);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->flags);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->time);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->tv);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->sdirk);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->gcn);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->nwt);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->psd);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->fluid);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->flow);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->source);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->buffer);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->pert);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->str);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->func);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->boost);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->dmda);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->recons);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->wbc);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->view);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->frc1d);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->frc2d);CHKERRQ(ierr);
    ierr = PetscMalloc1(1,&solve->frc3d);CHKERRQ(ierr);
    
    ierr = PetscStrcpy(solve->str->Input,args[1]);CHKERRQ(ierr);
    
    /*******************************************************************
     *                     Inicitalize Simulation
     ******************************************************************/
    
    ierr = read_parameters(solve);CHKERRQ(ierr);
    ierr = set_parameters(solve);CHKERRQ(ierr);
    ierr = set_fluid(solve);CHKERRQ(ierr);
    ierr = set_flow(solve);CHKERRQ(ierr);
    ierr = initialize(solve);CHKERRQ(ierr);
    
    func   = solve->func;
    dmda   = solve->dmda;
    flags  = solve->flags;
    //buffer = solve->buffer;
    //pert   = solve->pert;
    //view   = solve->view;
    
    /*******************************************************************
     *        Create Viewer to save solution during simulation
     ******************************************************************/
    
    ierr = create_viewers(solve);CHKERRQ(ierr);
    
    /*******************************************************************
     *                         Create DMDA
     ******************************************************************/
    
    ierr = (*func->dmda)(&dmda->DAInv,dmda->ne,3,solve);CHKERRQ(ierr);
    
    if ( dmda->_DIM_ == 2 ) {
		ierr = DMDASetFieldName(dmda->DAInv,0,"rho");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAInv,1,"rhou");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAInv,2,"rhov");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAInv,3,"rhoe");CHKERRQ(ierr);
	}else if ( dmda->_DIM_ == 3 ) {
		ierr = DMDASetFieldName(dmda->DAInv,0,"rho");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAInv,1,"rhou");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAInv,2,"rhov");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAInv,3,"rhow");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAInv,4,"rhoe");CHKERRQ(ierr);
	}
    
    ierr = (*func->dmda)(&dmda->DAVisc,dmda->ne,dmda->ghv,solve);CHKERRQ(ierr);
    
    ierr = (*func->dmda)(&dmda->DAFlux,dmda->ne,3,solve);CHKERRQ(ierr);
    
    ierr = (*func->dmda)(&dmda->DAOut,dmda->nout,1,solve);CHKERRQ(ierr);
    
    if ( dmda->_DIM_ == 2 ) {
		ierr = DMDASetFieldName(dmda->DAOut,0,"rho");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,1,"p");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,2,"u");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,3,"v");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,4,"T");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,5,"Ma");CHKERRQ(ierr);
	}else if ( dmda->_DIM_ == 3 ) {
		ierr = DMDASetFieldName(dmda->DAOut,0,"rho");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,1,"p");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,2,"u");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,3,"v");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,4,"w");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,5,"T");CHKERRQ(ierr);
		ierr = DMDASetFieldName(dmda->DAOut,6,"Ma");CHKERRQ(ierr);
	}
    
    /*******************************************************************
     *                      Create Vectos
     ******************************************************************/
    
    ierr = create_vectors(solve);CHKERRQ(ierr);
    
    
    /*******************************************************************
     *                    Call Initial Condition
     ******************************************************************/
    
    ierr = (*func->initialCondition)(solve);CHKERRQ(ierr);
    
    //ierr = total_variation(solve);CHKERRQ(ierr);
    
    /*******************************************************************
     *                       Set Source Terms
     ******************************************************************/
    
    //if( buffer->use ){
        //ierr = setBuffer(solve);CHKERRQ(ierr);
    //}
    
    //if( pert->use ){
        //ierr = computePerturbation(solve);CHKERRQ(ierr);
    //}
    
    /*******************************************************************
     *                      Time Integration 
     ******************************************************************/
    
    //ierr = PetscTime(&t1);CHKERRQ(ierr);
    
    //ierr = (*func->timeIntegration)(solve);CHKERRQ(ierr);
    
    output_solution_2d(NULL,0,0.0,dmda->qc,solve);
    
    //ierr = PetscTime(&t2);CHKERRQ(ierr);
    //ierr = PetscViewerASCIIPrintf( view->viewerCPUTime, "%.15g\n",t2-t1 );CHKERRQ(ierr);
    //ierr = saveSolution(1,dmda->qc,solve);CHKERRQ(ierr);

    //ierr = total_variation(solve);CHKERRQ(ierr);
    
    /*******************************************************************
     *                     Free Memory
     ******************************************************************/
    
    ierr = delete_viewers(solve);CHKERRQ(ierr);
    
    ierr = delete_vectors(solve);CHKERRQ(ierr);
    
    ierr = DMDestroy(&dmda->DAInv);CHKERRQ(ierr);
    ierr = DMDestroy(&dmda->DAVisc);CHKERRQ(ierr);
    ierr = DMDestroy(&dmda->DAFlux);CHKERRQ(ierr);
    ierr = DMDestroy(&dmda->DAOut);CHKERRQ(ierr);
    
    ierr = PetscFree(solve->prmt);CHKERRQ(ierr);
    ierr = PetscFree(solve->flags);CHKERRQ(ierr);
    ierr = PetscFree(solve->time);CHKERRQ(ierr);
    ierr = PetscFree(solve->tv);CHKERRQ(ierr);
    ierr = PetscFree(solve->sdirk);CHKERRQ(ierr);
    ierr = PetscFree(solve->gcn);CHKERRQ(ierr);
    ierr = PetscFree(solve->nwt);CHKERRQ(ierr);
    ierr = PetscFree(solve->psd);CHKERRQ(ierr);
    ierr = PetscFree(solve->fluid);CHKERRQ(ierr);
    ierr = PetscFree(solve->flow);CHKERRQ(ierr);
    ierr = PetscFree(solve->source);CHKERRQ(ierr);
    ierr = PetscFree(solve->buffer);CHKERRQ(ierr);
    ierr = PetscFree(solve->pert);CHKERRQ(ierr);
    ierr = PetscFree(solve->str);CHKERRQ(ierr);
    ierr = PetscFree(solve->func);CHKERRQ(ierr);
    ierr = PetscFree(solve->boost);CHKERRQ(ierr);
    ierr = PetscFree(solve->dmda);CHKERRQ(ierr);
    ierr = PetscFree(solve->recons);CHKERRQ(ierr);
    ierr = PetscFree(solve->wbc);CHKERRQ(ierr);
    ierr = PetscFree(solve->view);CHKERRQ(ierr);
    ierr = PetscFree(solve->frc1d);CHKERRQ(ierr);
    ierr = PetscFree(solve->frc2d);CHKERRQ(ierr);
    ierr = PetscFree(solve->frc3d);CHKERRQ(ierr);
    ierr = PetscFree(solve);CHKERRQ(ierr);
    
    ierr = PetscFinalize();
    return 0;
}
