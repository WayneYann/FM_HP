#include "FlameMaster.h"
#include "TCountDiffFlamePhys.h"
#include "TCountDiffPhysEigen.h"
//#include "TCountDiffFlameSim.h"
#include "TCountDiffFlameMix.h"
#include "TCountDiffFlameCont.h"
#include "TUnstrPremFlamePhys.h"
#include "TCountPremFlamePhys.h"
#include "TCountPremFlameSim.h"
#include "TFlameSheet.h"
#include "T0DIsoChor.h"
#include "TTransFlamelet.h"
//#include "TTrans1DIsoChor.h"
//#include "TTransCoalParticle.h"

int main( int argc, char *argv[] )
{

#ifdef __MWERKS__
	argc = ccommand(&argv);
#endif
	FirstInputPtr input = new FirstInput( argc, argv );
   	if ( !input ) FatalError( "memory allocation of FirstInputPtr failed" );
	//	input->PrintAdditionalData();
	if ( input->fFlameType == kStartUpDiffusion ) {
		TFlameSheetPtr flame = new TFlameSheet( input );
 	  	if ( !flame ) FatalError( "memory allocation of TFlameSheetPtr failed" );
		flame->GetInputData()->PrintAdditionalData();
		flame->GetSolver()->Solve( flame );
		delete flame;
	}
	else if ( input->fFlameType == kDiffusionPhys )  {
		TCountDiffFlamePhysPtr flame = new TCountDiffFlamePhys( input );
  	 	if ( !flame ) FatalError( "memory allocation of TDiffusionFlamePtr failed" );
		flame->GetSolver()->Solve( flame );
		delete flame;
	}
	else if ( input->fFlameType == kDiffPhysEigen )  {
		TCountDiffPhysEigenPtr flame = new TCountDiffPhysEigen( input );
  	 	if ( !flame ) FatalError( "memory allocation of TCountDiffPhysEigenPtr failed" );
		flame->GetSolver()->Solve( flame );
		delete flame;
	}
	else if ( input->fFlameType == kUnstrPremPhys )  {
		TUnstrPremFlamePhysPtr flame = new TUnstrPremFlamePhys( input );
  	 	if ( !flame ) FatalError( "memory allocation of TPremixedFlamePtr failed" );
		flame->GetSolver()->Solve( flame );
		delete flame;
	}
	else if ( input->fFlameType == kCountPremPhys )  {
		TCountPremFlamePhysPtr flame = new TCountPremFlamePhys( input );
  	 	if ( !flame ) FatalError( "memory allocation of TDiffusionFlamePtr failed" );
		flame->GetSolver()->Solve( flame );
		delete flame;
	}
	else if ( input->fFlameType == kCountPremSim )  {
		TCountPremFlameSimPtr flame = new TCountPremFlameSim( input );
  	 	if ( !flame ) FatalError( "memory allocation of TDiffusionFlamePtr failed" );
		flame->GetSolver()->Solve( flame );
		delete flame;
	}
	else if ( input->fFlameType == kStartUpDiffusion2 ) {
		TCountDiffFlameSimPtr flame = new TCountDiffFlameSim( input );
   		if ( !flame ) FatalError( "memory allocation of TCountDiffFlameSim failed" );
		flame->GetSolver()->Solve( flame );
		delete flame;
	}
/* 	else if ( input->fFlameType == kCountDiffCont ) { */
/* //		input->fFlameType = kStartUpDiffusion2; */
/* 		input->fVariablesWithoutSpecies = 4; */
/* 		TCountDiffFlameContPrePtr flamePre = new TCountDiffFlameContPre( input ); */
/*    		if ( !flamePre ) FatalError( "memory allocation of TCountDiffFlameContPre failed" ); */
/* 		flamePre->GetSolver()->Solve( flamePre ); */
/* 		Double		dTds = flamePre->GetdTds(); */
/* 		Double		dAInvds = flamePre->GetdAInvds(); */
/* 		Double		strainRate = flamePre->GetStrainRate(); */
/* 		TSolutionPtr	sol = flamePre->GetSolution(); */
/* 		 */
/* 		delete flamePre; */
/*  */
/* //		input->fFlameType = kCountDiffCont; */
/* 		input->fVariablesWithoutSpecies = 5; */
/* 		TCountDiffFlameContPtr flame = new TCountDiffFlameCont( input, sol, strainRate ); */
/* 		delete sol; */
/* 		flame->SetdTds( dTds ); */
/* 		flame->SetdAInvds( dAInvds ); */
/*    		if ( !flame ) FatalError( "memory allocation of TCountDiffFlameCont failed" ); */
/* 		flame->GetSolver()->Solve( flame ); */
/* 		delete flame; */
/* 	} */
	else if ( input->fFlameType == kCountDiffMix ) {
		TCountDiffFlameMixPtr flame = new TCountDiffFlameMix( input );
   		if ( !flame ) FatalError( "memory allocation of TCountDiffFlameMix failed" );
		flame->GetSolver()->Solve( flame );
		delete flame;
	}
	else if ( input->fFlameType == kHomoIsoChor || input->fFlameType == kHomoIsoBar ) {
		T0DIsoChorPtr flame = new T0DIsoChor( input );
/*		gFlame = flame;
		if ( atexit( ExitFunc ) ) {
			cerr << "#warning: installation of ExitFunc failed" << NEWL;
		}*/
   		if ( !flame ) FatalError( "memory allocation of T0DIsoChor failed" );
		flame->Solve();
		delete flame;
	}
	else if ( input->fFlameType == kTransFlamelet ) {
		TTransFlameletPtr flame = new TTransFlamelet( input );
   		if ( !flame ) FatalError( "memory allocation of TTransFlamelet failed" );
		flame->Solve();
		delete flame;
	}
/* 	else if ( input->fFlameType == kTrans1DIsoChor ) { */
/* 		TTrans1DIsoChorPtr flame = new TTrans1DIsoChor( input ); */
/*    		if ( !flame ) FatalError( "memory allocation of TTrans1DIsoChor failed" ); */
/* 		flame->Solve(); */
/* 		delete flame; */
/* 	} */
/* 	else if ( input->fFlameType == kTransCoalParticle ) { */
/* 		TTransCoalParticlePtr flame = new TTransCoalParticle( input ); */
/*    		if ( !flame ) FatalError( "memory allocation of TTransCoalParticle failed" ); */
/* 		flame->Solve(); */
/* 		delete flame; */
/* 	} */
	else { 
		cerr << "error: wrong flame type specified: " << input->fFlameType << NEWL;
		exit(2);
	}
	delete input;
	fputc( '\a', stderr );
	return 0;
}

