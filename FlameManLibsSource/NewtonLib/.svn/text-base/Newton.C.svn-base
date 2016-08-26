#include "Newton.h"

#undef DEBUGMODI

#undef HIGHSOPH
#define NEWSCALEDSOLUTION

#define TIMEANDDAMP

AMPrintOptionsPtr gPrnt;


TBVPSolver::TBVPSolver( TBVPInputPtr input )
{
  bt = new TNewton( input );
  if ( !bt ) NewtonFatalError( "memory allocation of TNewton failed" );
  
  // Take into account the possibility to have Time dependancy 
  // and Damping after switching off time dependance

  if ( bt->GetDampFlag() && !bt->GetTimedepFlag()){
    damp = new TDamp( bt, input );
    if ( !damp ) NewtonFatalError( "memory allocation of TDamp failed" );
    time = NULL;
  }   

  else if ( bt->GetTimedepFlag() && !bt->GetDampFlag()){  
    time = new TTime( bt, input );		
    if ( !time ) NewtonFatalError( "memory allocation of TTime failed" );
    damp = NULL;
  }   
  
  else if ( bt->GetTimedepFlag() && bt->GetDampFlag()){    
    time = new TTime( bt, input );		
    if ( !time ) NewtonFatalError( "memory allocation of TTime failed" );
    damp = new TDamp( bt, input );
    if ( !damp ) NewtonFatalError( "memory allocation of TDamp failed" );
    bt->SetDampFlag( FALSE ); //no damping during time dependant computations
  }   

  else{
    time = NULL;
    damp = NULL;
  }

  if ( bt->GetContinuationFlag() ) {
    contin = new TContinuation( bt, input );	
    if ( !contin ) NewtonFatalError( "memory allocation of TContinuation failed" );
  }
  else {
    contin = NULL;
  }
}

TBVPSolver::~TBVPSolver( void )
{
//	delete problem;
	delete contin;
	delete time;
	delete damp;
	delete bt;
}

int TBVPSolver::Solve( void *object )
{
	int         	err = 0;                //  error handling
	TGridPtr		currentGrid = bt->GetGrid()->GetCurrentGrid();

	AMPrintOptions prnt;
	gPrnt = &prnt;
	DefaultAMPOpts( gPrnt );
	FILE 	*fpb, *fpa;
	FILE	*fp;
	int		counter;
#if defined (applec) || defined (powerc)
	Double value;			
#else
	Double	value;
#endif

	Test_convergence( this, object );

//  begin computation
	do {            //  continuation loop
		if ( bt->GetContinuationFlag() ){
			contin->Contin_step( this, damp );
		}
		if ( time ) {
			time->InitTimeMats( bt, object );
//			time->InitYPrime( bt, object );
		}
		do {        //  Newton loop
			bt->IncNIter();
	
			bt->CheckGridding( time );

//  Newton step
			err = bt->Newton_step( damp, time, object );
			if ( err ) break;

			if ( bt->IsWriteEverySolution() ) {
				counter = ( int ) ( modf( bt->GetNIter()/10.0, &value ) * 10.0 + 0.1 );
				sprintf( bt->GetOutFileBuff(), "%sdat%d.dout", bt->GetOutputPath(), counter );
				if ( ( fp = fopen( bt->GetOutFileBuff(), "w") ) == NULL ) {
					cerr << "unable to open file " << bt->GetOutFileBuff() << NEWL;
					exit( 2 );
				}
				bt->PrintSolution( currentGrid->GetX()->vec, currentGrid->GetY()->mat, bt->GetVariableNames( object ), fp );
				fclose( fp );
				if ( counter > 10 ) counter = 1;
			}
			
			Test_convergence( this, object );

			Monitor( bt, contin, damp, time );

			if ( bt->IsGriddingStep() && bt->IsWatchGridding() ) {
				sprintf( bt->GetOutFileBuff(), "%soldGrid.dout", bt->GetOutputPath() );
				fpb = fopen( bt->GetOutFileBuff(), "w" );
				bt->PrintSolution( currentGrid->GetX()->vec, currentGrid->GetY()->mat, bt->GetVariableNames( object ), fpb );
				fclose( fpb );
				if ( time ) {
//					sprintf( bt->GetOutFileBuff(), "%yPrimeV.dout", bt->GetOutputPath() );
//					fpb = fopen( bt->GetOutFileBuff(), "w" );
//					bt->PrintSolution( currentGrid->GetX()->vec, time->GetYPrime()->mat, bt->GetVariableNames( object ), fpb );
//					fclose( fpb );
				}
			}
			if ( !bt->GetLeaveNewton() ) {
				bt->GetGrid()->new_grid( this, object );
			}
			if ( bt->IsGriddingStep() && bt->IsWatchGridding()  ) {
				bt->GetUTFuncs()->WriteOutput( object, NULL, "grid" );
				sprintf( bt->GetOutFileBuff(), "%snewGrid.dout", bt->GetOutputPath() );
				fpa = fopen( bt->GetOutFileBuff(), "w" );
				bt->PrintSolution( currentGrid->GetX()->vec, currentGrid->GetY()->mat, bt->GetVariableNames( object ), fpa );
				fclose( fpa );
				if ( time ) {
//					sprintf( bt->GetOutFileBuff(), "%yPrimeN.dout", bt->GetOutputPath() );
//					fpa = fopen( bt->GetOutFileBuff(), "w" );
//					bt->PrintSolution( currentGrid->GetX()->vec, time->GetYPrime()->mat, bt->GetVariableNames( object ), fpa );
//					fclose( fpa );
				}
			}
            
			if ( !bt->IsGriddingStep() && bt->GetTimedepFlag() && time->GetAdjustTime() 
				&& bt->GetNIter() > 1 ) {

			  // Check that we do not exit on a modified Newton
			  if (!(bt->IsFullStep ())) {
			    bt->ForceFullStep ();
			  }
			  else {
			    if ( bt->GetLeaveNewton() ) {  // this time step converged
			      ///  fprintf( stderr, "unset leaveNewton1\n" );
			      bt->UnSetLeaveNewton();
			      time->UnSetTimeConverged();
			      time->AdjustTimeStep( this, FALSE, object ); // prepare for new time step
			    }
			    else {
			      time->AdjustTimeStep( this, TRUE, object ); // prepare for new time step
			    }

			  }
			}
#ifdef HP
			if ( gExit ) {
				DoExit( bt, object );
			}
#elif defined SUNGXX
			if ( gExit ) {
				DoExit( bt, object );
			}
#endif
		} while ( !bt->GetLeaveNewton() );

		if ( bt->GetConvergeNewton() ) {
		  // Check that we do not exit on an modified Newton
		  if (!(bt->IsFullStep ())) {
		    bt->ForceFullStep ();
		    bt->UnSetLeaveContin ();
		  }
		  else {
			cerr << "* newton iteration converged *" << NEWL << NEWL;
			bt->PostConvergence( object );
		  }
		}
		else {
			cerr << "# no convergence" << NEWL;
			bt->WriteOutput( object, NULL, "noC" );
			bt->PostConvergence( object );
		}
	} while ( !bt->GetLeaveContin() );
	cerr << "done" << NEWL;
	return 0;
}

void TBVPSolver::UpdateAllDimensions( int nOfPoints )
{
	bt->AdjustNewtonDimensions( nOfPoints );
	bt->GetGrid()->UpdateGridDimension();
	
	if ( bt->GetDampFlag() ) {
		damp->AdjustDampDimension( nOfPoints );
	}
	if ( bt->GetTimedepFlag() ) {
		time->AdjustTimeDimension( nOfPoints );
	}
	if ( bt->GetContinuationFlag() ) {
		contin->AdjustContinDimension( nOfPoints );
	}
}

void TNewton::InitNewton( TBVPInputPtr input )
{
  max_gridpoints = input->fMaxGridPoints;
  n_equations = input->fNVariables;

  if ( !( max_gridpoints % 2 ) ) ++max_gridpoints;

  initialGridPoints = ( int ) MIN( input->fInitialGridPoints, max_gridpoints );
  //	initialGridPoints = ( int ) MAX( initialGridPoints, 5 );

  timedep_flag = input->fTimeDepFlag;
  
  // Damp_flag is set in any case and set to FALSE during 
  // time-dependant computations 
  damp_flag = input->fDampFlag;
  
// #ifdef TIMEANDDAMP
// 	damp_flag = input->fDampFlag;
// #else
// 	if ( !timedep_flag ) {
// 		damp_flag = input->fDampFlag;
// 	}
// 	else {
// 		if ( input->fDampFlag ) {
// 			fprintf( stderr, "switch off damping\n" );
// 			damp_flag = FALSE;
// 		}
// 	}
// #endif

	continuation_flag = input->fContinFlag;
	fWithModified = input->fUseModifiedNewton;
	fUseNumericalJac = input->fUseNumericalJac;
	fUseSecOrdJac = input->fUseSecOrdJac;
	fModCounter = 0;
	fWriteFullRes = input->fWriteFullRes;
	fWriteBT = input->fWriteBT;
	fWriteResiduum = input->fWriteResiduum;
	fWatchGridding = input->fWatchGridding;

	fWriteEverySolution = input->fWriteEverySolution;
	fOutputPath = GetFullPath( input->fOutputPath, kPath );
//	cerr << "outpath is '" << fOutputPath << "'" << NEWL;

	fOutFile = new char[strlen( fOutputPath ) + 128];

	delta_new_grid = input->fDeltaNewGrid;
	parameter = NULL;

	TolRes = ( input->fTolRes > 0.0 ) ? input->fTolRes : 1.0e-15;
	TolDy = ( input->fTolDy > 0.0 ) ? input->fTolDy : 1.0e-5;

    max_iter = input->fMaxIter;
//	left = input->fLeft;
//	right = input->fRight;

	dy = NewMatrix( n_equations, max_gridpoints, kColumnPointers );
	dySaved = NewMatrix( n_equations, max_gridpoints, kColumnPointers );
    fA = NewTensor( max_gridpoints, n_equations, n_equations, kColumnPointers );
    fB = NewTensor( max_gridpoints, n_equations, n_equations, kColumnPointers );
    fC = NewTensor( max_gridpoints, n_equations, n_equations, kColumnPointers );
    ip = New1DIntArray( max_gridpoints * n_equations );
	fJacobiScaler = NewMatrix( n_equations, max_gridpoints, kColumnPointers );
	fYMax = NewVector( n_equations );
	fResiduals = NewVector( n_equations );
	fNodeInfo = new NodeInfo();
	if ( !fNodeInfo ) NewtonFatalError( "memory allocation of NodeInfo failed" );
    grid = new TAdaptiveGrid( this, input );
	if ( !grid ) NewtonFatalError( "memory allocation of TAdaptiveGrid failed" );
	dx = ( GetRight() - GetLeft() ) / ( initialGridPoints + 1 );
	
//	initialization	
	if ( input->fInitialEquations > 0 && input->fInitialEquations < n_equations ) {
		SetNOfEquations( input->fInitialEquations );
		cerr << "start to solve " << input->fInitialEquations << " equations" << NEWL;
	}
    n_iter = 0;
    converg_newton = 0;
    leave_newton = 0;
    leave_contin = 0;
	fModified = FALSE;
    normOldRes = 10.0;
    normRes = 10.0;
    normOldDy = 10.0;
    normDy = 10.0;
	InitModifiedConstants();
	fUtFuncs = new UTFuncs();
	
#	ifdef HP
	InstallInterruptHP();
#	endif
}

TNewton::~TNewton( void )
{
    delete fUtFuncs;
    delete grid;
	delete fNodeInfo;
	DisposeVector( fResiduals );
	DisposeVector( fYMax );
	DisposeMatrix( fJacobiScaler );
	Free1DIntArray( ip );
    DisposeTensor( fC );
    DisposeTensor( fB );
    DisposeTensor( fA );
    DisposeMatrix( dySaved );
    DisposeMatrix( dy );
	delete fOutFile;
	delete fOutputPath;
} 

int TNewton::Newton_step( TDampPtr damp, TTimePtr time, void *object )
{
        int 		err = 0;
        TGridPtr	currentGrid = grid->GetCurrentGrid();
        AMPrintOptions  prnt;

	FILE	*fp;
    
	if ( !forcefullstep && CheckModified() ) {
		ClearMatrix( dy );
		BackSolve( dy, object );
/*		if ( damp_flag && grid->IsFine() ) {*/
/*			damp->Update_Lambda( this, time, object );*/
/*			if ( IsWriteBT() ) {*/
/*				fprintf( fp, "\nLambda = %g\n", damp->GetLambda() );*/
/*				fprintf( stderr, "Lambda out\n" );*/
/*			}*/
/*		}*/


		copy_mat( dySaved, currentGrid->GetY() ); // save solution in case modified step is bad
		if ( UpdateSolution( object, damp, time ) ) {
			fullstep = TRUE;
		}
		else {
			grid->SetSolutionScaler();
			fullstep = !AcceptModified();
		}
		if ( fullstep ) {
			copy_mat( currentGrid->GetY(), dySaved ); // restore solution
			fUtFuncs->PostIter( object );
		}
		else {
			--n_iter;
		}
	}
	else {
		fullstep = TRUE;
		fModified = FALSE;
	}

	if ( fullstep ) {
	    ClearMatrix( dy );
	    ClearTensor( fA );
	    ClearTensor( fB );
	    ClearTensor( fC );

	    if ( fUseNumericalJac ) 
	      {
		UpdateJacRHSNum( object ); 
	      }
	    else 
	      {
		UpdateJacRHS( object );
	      }

//		if ( timedep_flag ) Update_timedepJac( time );

		
#ifdef NEWSCALEDSOLUTION
		ScaleSystem();
#endif

		if ( IsWriteBT() && grid->IsFine() ) {
			DefaultAMPOpts( &prnt );
			sprintf( GetOutFileBuff(), "%sbt.tout", GetOutputPath() );
			if ( !( fp = fopen( GetOutFileBuff(), "w" ) ) ) {
				fprintf( stderr, "#error: unable to open output file '%s'\n"
							, GetOutFileBuff() );
				exit(2);
			}
			fprintf( fp, "dx = %g", dx );

			prnt.format = "%12.4E";

			prnt.title = "\nTensor a, GridPoint no. 0";
			fprintf( fp, "\nnIter = %d", n_iter );
			PrintMatrix( SubMatrix( 0, fA ), &prnt, fp );
			fprintf( stderr, "Tensor a out\n" );
			
			prnt.title = "\nTensor b, GridPoint no. 0";
			fprintf( fp, "\nnIter = %d", n_iter );
			PrintMatrix( SubMatrix( 0, fB ), &prnt, fp );
			fprintf( stderr, "Tensor b out\n" );

			prnt.title = "\nTensor c, GridPoint no. 0";
			fprintf( fp, "\nnIter = %d", n_iter );
			PrintMatrix( SubMatrix( 0, fC ), &prnt, fp );
			fprintf( stderr, "Tensor b out\n" );

			prnt.title = "\nMatrix RHS";
			fprintf( fp, "\nnIter = %d", n_iter );
			PrintMatrix( dy, &prnt, fp );
			fprintf( stderr, "RHS out\n" );

//			prnt.title = "\nMatrix JacobiScaler";
//			fprintf( fp, "\nnIter = %d", n_iter );
//			PrintMatrix( fJacobiScaler, &prnt, fp );

			fflush( fp );
		}

		err = decbt( fA, fB, fC, ip );
		if ( err ) {
			fprintf( stderr, "# error in decbt:" ); 
			if ( err == -1 ) { fprintf( stderr, "input value of order of blocks or number of blocks illegal\n" ); }
			else if ( err == -2 ) { fprintf( stderr, "invalid partition scheme (kColumnPointers needed)\n" ); }
			else { fprintf(stderr, "singular matrix found in %dth diagonal block\n", err ); }
//			exit( 2 );
//			UnSetConvergeNewton();
//			SetLeaveNewton();
//			return err;
		}

		err = solbt( fA, fB, fC, dy, ip );
		if ( err ) { fprintf(stderr, "# error in solbt\n"); exit( 2 ); }

		if ( damp_flag && grid->IsFine() ) {
			damp->Update_Lambda( this, time, object );
			if ( IsWriteBT() ) {
				fprintf( fp, "\nLambda = %g\n", damp->GetLambda() );
				fprintf( stderr, "Lambda out\n" );
			}
		}
	
		if ( IsWriteBT() ) {
			prnt.title = "\nMatrix yBeforeUpdate";
			fprintf( fp, "\nnIter = %d", n_iter );
			PrintMatrix( grid->GetCurrentGrid()->GetY(), &prnt, fp );
			fprintf( stderr, "Matrix yBeforeUpdate out\n" );
		}
	
		if ( UpdateSolution( object, damp, time ) ) {
			UnSetConvergeNewton();
			SetLeaveNewton();
			return 1;
		}
	//	time->InitYPrime( this, object );
	/*	if ( timedep_flag ) {*/
	/*		time->UpdateYPrime( this );*/
	/*	}*/
	
		grid->SetSolutionScaler();
#ifdef DEBUGMODI
		FILE        *fpy;
		char        fgame[32];
	
		sprintf( fgame, "scalerNew%d.tout", GetNIter() );
		fpy = fopen( fgame, "w" );
		gPrnt->title = "\nVector scaler";
		PrintVector( grid->GetSolutionScaler(), gPrnt, fpy );
		fclose(fpy);
#endif
	}

	if ( fWithModified ) {
		fBMod = ScaledMaxNorm( dy, grid->GetSolutionScaler() );
	}

	if ( IsWriteBT() ) {
		fprintf( fp, "\nnIter = %d", n_iter );
		prnt.title = "\nMatrix dy";
		PrintMatrix( dy, &prnt, fp );
		fprintf( stderr, "Matrix dy out\n" );

		fprintf( fp, "\nnIter = %d", n_iter );
		prnt.title = "\nMatrix y";
		PrintMatrix( grid->GetCurrentGrid()->GetY(), &prnt, fp );
		fprintf( stderr, "Matrix y out\n" );

		fclose( fp );

		fprintf( stderr, "bt dumped for nIter = %d\n", n_iter );
//		exit(2);
	}

	// Reset forcefullstep to FALSE
	forcefullstep = FALSE;

	return 0;
}

void TNewton::InitModifiedConstants( void )
{
//	following Smooke's approach fAMod has to be equal to 0.5
//	for a higher rate of convergence , I set fAMod equal to 1.0
	fAMod = 1.0;
	fCMod = 1.0;
	fBMod = 0.0;
	
	fModified = FALSE;
}

Flag TNewton::AcceptModified( void )
{

	if ( fBMod == 0.0 ) return FALSE;

#ifdef HIGHSOPH
	static Double	h = 0.5;
	Double			locAMod = fAMod;
	if ( ScaledMaxNorm( dy, grid->GetSolutionScaler() ) < h * fBMod * fAMod ) {
		fAMod = h * locAMod * ( fCMod + 0.5 * h * locAMod );	// means fCModOld
		fCMod += h * locAMod;
		return TRUE;
	}
	else {
		InitModifiedConstants();
		return FALSE;
	}
#else
#ifdef DEBUGMODI
    FILE        *fpy;
    char        fgame[32];

    sprintf( fgame, "scalermod%d.tout", GetNIter() );
    fpy = fopen( fgame, "w" );
    gPrnt->title = "\nVector scaler";
    PrintVector( grid->GetSolutionScaler(), gPrnt, fpy );
    fclose(fpy);
#endif
	if ( ScaledNorm( dy, grid->GetSolutionScaler() ) < normDy ) {
//		fprintf( stderr, "normOldDy = %g\tScaledNorm = %g\n", normDy, ScaledNorm( dy, grid->GetSolutionScaler() ) );
		return TRUE;
	}
	else {
		fModified = FALSE;
		return FALSE;
	}
#endif
}

int TNewton::UpdateSolution( void *object, TDampPtr damp, TTimePtr time )
{
	Update_Y( damp, time );
	return fUtFuncs->PostIter( object );
}

int TNewton::UpdateJacRHS( void *object )
{
	TGridPtr			currentGrid = grid->GetCurrentGrid();
	int             	k;
	int 				N = currentGrid->GetNGridPoints();
	enum            	{ first, rest, last };

#if defined (applec) || defined (powerc)
	RotateCursor( 32 * n_iter );
#endif

	// first point
	SetNodeInfo( object, 0 );
	fUtFuncs->JacobiFirstPoint( object, fNodeInfo );
	fUtFuncs->RHSFirstPoint( object, fNodeInfo, kDoNothing );

// all points except boundaries
	for ( k = 1; k < N-1; ++k ){
#if defined (applec) || defined (powerc)
	RotateCursor( 32 );
#endif
	//cerr << "set jac and rhs for point no. " << k << NEWL;
		SetNodeInfo( object, k );
		fUtFuncs->JacobiRestPoint( object, fNodeInfo );
		fUtFuncs->RHSRestPoint( object, fNodeInfo, kDoNothing );
	}

// last point
	SetNodeInfo( object, N-1 );
	fUtFuncs->JacobiLastPoint( object, fNodeInfo );
	fUtFuncs->RHSLastPoint( object, fNodeInfo, kDoNothing );

    return 0;
}

int TNewton::UpdateJacRHSNum( void *object )
{
  int i, j, k;
  TGridPtr currentGrid = grid->GetCurrentGrid ();
  double **Y = currentGrid->GetY ()->mat;
  double *locY;
  double savedY;
  double deltaY;
  int N = currentGrid->GetNGridPoints ();
  int nVariables = GetNVariables ();
  int nEquations = GetNEquations ();
  enum { first, rest, last };
  double alpha = 1.0e-6;
  double beta = 1.0e-14;

  UpdateRHS (dySaved, object);

  // First point (Left)
  locY = Y[0];

  for (j = 0; j < nVariables; ++j) {
    deltaY = alpha * locY[j] + beta;
    savedY = locY[j];

    // Upwind First Order 
    locY[j] += deltaY;

    fUtFuncs->UpdateLeftBoundary (object);
    SetNodeInfo (object, 0);
    fUtFuncs->RHSFirstPoint (object, fNodeInfo, kUpdate);
    for (i = 0; i < nEquations; ++i) {
      fNodeInfo->a[j][i] -=
	(fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
    }

    SetNodeInfo (object, 1);
    fUtFuncs->RHSFirstPoint (object, fNodeInfo, kDoNothing);
    for (i = 0; i < nEquations; ++i) {
      fNodeInfo->c[j][i] -=
	(fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
    }

    locY[j] = savedY;

    // Centered Second Order
    if (fUseSecOrdJac) {
      locY[j] -= deltaY;
      if (locY[j] > 0) {

	fUtFuncs->UpdateLeftBoundary (object);
	SetNodeInfo (object, 0);
	fUtFuncs->RHSFirstPoint (object, fNodeInfo, kUpdate);
	for (i = 0; i < nEquations; ++i) {
	  fNodeInfo->a[j][i] -=
	    (fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
	  fNodeInfo->a[j][i] /= 2.0;
	}

	SetNodeInfo (object, 1);
	fUtFuncs->RHSFirstPoint (object, fNodeInfo, kDoNothing);
	for (i = 0; i < nEquations; ++i) {
	  fNodeInfo->c[j][i] -=
	    (fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
	  fNodeInfo->c[j][i] /= 2.0;
	}
      }
      locY[j] = savedY;
    }
  }

  fUtFuncs->UpdateLeftBoundary (object);
  SetNodeInfo (object, 0);
  fUtFuncs->RHSFirstPoint (object, fNodeInfo, kUpdate);

  // all points except boundaries
  for (k = 1; k < N - 1; ++k) {
    locY = Y[k];

    for (j = 0; j < nVariables; ++j) {

      deltaY = alpha * locY[j] + beta;
      savedY = locY[j];

      // Upwind First Order
      locY[j] += deltaY;

      SetNodeInfo (object, k);
      fUtFuncs->RHSRestPoint (object, fNodeInfo, kUpdate);
      for (i = 0; i < nEquations; ++i) {
	fNodeInfo->a[j][i] -=
	  (fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
      }

      SetNodeInfo (object, k - 1);
      fUtFuncs->RHSRestPoint (object, fNodeInfo, kDoNothing);
      for (i = 0; i < nEquations; ++i) {
	fNodeInfo->b[j][i] -=
	  (fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
      }

      SetNodeInfo (object, k + 1);
      fUtFuncs->RHSRestPoint (object, fNodeInfo, kDoNothing);
      for (i = 0; i < nEquations; ++i) {
	fNodeInfo->c[j][i] -=
	  (fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
      }
      locY[j] = savedY;

      // Centered Second Order
      if (fUseSecOrdJac) {
	locY[j] -= deltaY;
	if (locY[j] > 0) {

	  SetNodeInfo (object, k);
	  fUtFuncs->RHSRestPoint (object, fNodeInfo, kUpdate);
	  for (i = 0; i < nEquations; ++i) {
	    fNodeInfo->a[j][i] -=
	      (fNodeInfo->rhs[i] -
	       fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
	    fNodeInfo->a[j][i] /= 2.0;
	  }

	  SetNodeInfo (object, k - 1);
	  fUtFuncs->RHSRestPoint (object, fNodeInfo, kDoNothing);
	  for (i = 0; i < nEquations; ++i) {
	    fNodeInfo->b[j][i] -=
	      (fNodeInfo->rhs[i] -
	       fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
	    fNodeInfo->b[j][i] /= 2.0;
	  }

	  SetNodeInfo (object, k + 1);
	  fUtFuncs->RHSRestPoint (object, fNodeInfo, kDoNothing);
	  for (i = 0; i < nEquations; ++i) {
	    fNodeInfo->c[j][i] -=
	      (fNodeInfo->rhs[i] -
	       fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
	    fNodeInfo->c[j][i] /= 2.0;
	  }
	}
	locY[j] = savedY;
      }
    }
    SetNodeInfo (object, k);
    fUtFuncs->RHSRestPoint (object, fNodeInfo, kUpdate);

  }

  // Last point (Right)
  locY = Y[N - 1];
  for (j = 0; j < nVariables; ++j) {

    deltaY = alpha * locY[j] + beta;
    savedY = locY[j];

    // Upwind First Order
    locY[j] += deltaY;

    fUtFuncs->UpdateRightBoundary (object);
    SetNodeInfo (object, N - 1);
    fUtFuncs->RHSLastPoint (object, fNodeInfo, kUpdate);
    for (i = 0; i < nEquations; ++i) {
      fNodeInfo->a[j][i] -=
	(fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
    }

    SetNodeInfo (object, N - 2);
    fUtFuncs->RHSLastPoint (object, fNodeInfo, kDoNothing);
    for (i = 0; i < nEquations; ++i) {
      fNodeInfo->b[j][i] -=
	(fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
    }

    locY[j] = savedY;

    // Centered Second Order
    if (fUseSecOrdJac) {
      locY[j] -= deltaY;
      if (locY[j] > 0) {

	fUtFuncs->UpdateRightBoundary (object);
	SetNodeInfo (object, N - 1);
	fUtFuncs->RHSLastPoint (object, fNodeInfo, kUpdate);
	for (i = 0; i < nEquations; ++i) {
	  fNodeInfo->a[j][i] -=
	    (fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
	  fNodeInfo->a[j][i] /= 2.0;
	}

	SetNodeInfo (object, N - 2);
	fUtFuncs->RHSLastPoint (object, fNodeInfo, kDoNothing);
	for (i = 0; i < nEquations; ++i) {
	  fNodeInfo->b[j][i] -=
	    (fNodeInfo->rhs[i] - fNodeInfo->rhsSaved[i]) / (locY[j] - savedY);
	  fNodeInfo->b[j][i] /= 2.0;
	}
      }
      locY[j] = savedY;
    }
  }
  fUtFuncs->UpdateRightBoundary (object);
  SetNodeInfo (object, N - 1);
  fUtFuncs->RHSLastPoint (object, fNodeInfo, kUpdate);

  // restore rhs
  CopyMatrix (dySaved, dy);

  return 0;
}

int TNewton::UpdateRHS( MatrixPtr locDyMat, void *object )
{
    TGridPtr			currentGrid = grid->GetCurrentGrid();
    int         		k;
    int         		N = GetCurrentGridPoints();
    enum        		{ first, rest, last };
	Double				**locDy = locDyMat->mat;
	
#if defined (applec) || defined (powerc)
    RotateCursor( 32 * n_iter );
#endif

	ClearMatrix( locDyMat );
	
// first point
	SetNodeInfo( object, 0 );
	fNodeInfo->rhs = locDy[0];	// could be an other address than bt->dy
	fUtFuncs->RHSFirstPoint( object, fNodeInfo, kDoNothing );

// all points except boundaries
	for ( k = 1; k < N-1; ++k ){
#if defined (applec) || defined (powerc)
	RotateCursor( 32 );
#endif
		SetNodeInfo( object, k );
		fNodeInfo->rhs = locDy[k];	// could be an other address than bt->dy
		fUtFuncs->RHSRestPoint( object, fNodeInfo, kDoNothing );
	}

// last point
	SetNodeInfo( object, N-1 );
	fNodeInfo->rhs = locDy[N-1];	// could be an other address than bt->dy
	fUtFuncs->RHSLastPoint( object, fNodeInfo, kDoNothing );

//	if ( timedep_flag )  Update_timedepRHS( time );

    return 0;
}

int TNewton::CheckModified( void )
{
	fModified = TRUE;

	if ( !fWithModified ) {
		fModified = FALSE;
	}

/*	if ( damp_flag ) {
		fModified = FALSE;
	}*/

#ifndef HIGHSOPH
//	if ( GetTimedepFlag() ) {
//		if ( !( ( ( normRes / normOldRes < 1.0 && normRes < 1.0 ) || 
//			   ( normDy / normOldDy < 1.0 && normDy < 1.0 ) ) && fModCounter < 10 && n_iter > 1 ) ) {
//			fModified = FALSE;
//		}
//	}
//	else {
		if ( !( ( ( normRes / MAX(1e-30,normOldRes) < 0.5 && normRes < 1.0 ) || 
			   ( normDy / MAX(1e-30,normOldDy) < 0.5 && normDy < 1.0 ) ) && n_iter > 1 ) ) {
			fModified = FALSE;
		}
//	}
#endif

	if ( IsGriddingStep() ) {
		fModified = FALSE;
	}

	if ( n_iter - grid->GetIterOfGeneration() == 1 ) {
		fModified = FALSE;
	}

#ifndef HIGHSOPH
	if ( GetTimedepFlag() ) {
		if ( fModified ) {
//			++fModCounter;
		}
		else {
//			fModCounter = 0;
		}
	}
#endif
	return fModified;
}

/*void TNewton::Update_timedepJac( TTimePtr time )*/
/*{*/
/*	TGridPtr	locGrid = grid->GetCurrentGrid();	*/
/*	int 		i, k;*/
/*	int			nGridPoints = locGrid->GetNGridPoints();*/
/*	Double		deltaT = time->GetDeltaT();*/
/*	Double		deltaTMax = time->GetDeltaTMax();*/
/*	Double		hnennOverDeltaT;*/
/*	Double 		h;*/
/*	Double 		hm;*/
/*	Double 		hnenn;*/
/*	Double 		*x = locGrid->GetX()->vec;*/
/*	Double		**aMat = NULL;*/
/*	Double		***a = fA->tensor;*/
/*	Double		**locDy = dy->mat;*/
/*	Double		**yPrime = time->GetYPrime()->mat;*/
/*	*/
/*	if ( deltaT >= deltaTMax ) {*/
/*		if ( !time->GetTimeConverged() ) {*/
/*			cerr << "switch off time dependence" << NEWL;*/
/*			grid->UnSetSolHasNewGrid();*/
/*			time->SetTimeConverged();*/
/*		}*/
/*		return;*/
/*	}*/
/*	*/
/*    for ( k = 0; k < nGridPoints; ++k ){*/
/*		aMat = a[k];*/
/*	*/
/*        if ( k == 0 ){*/
/*            hm = x[0] - GetLeft();*/
/*            h = x[1]- x[0];*/
/*        }*/
/*        else if ( k == nGridPoints-1 ){*/
/*            hm = x[k] - x[k-1];*/
/*            h = GetRight() - x[k];*/
/*        }*/
/*        else{*/
/*            hm = x[k] - x[k-1];*/
/*            h = x[k+1] - x[k];*/
/*        }*/
/*		hnenn = h * hm * ( h + hm );*/
/* 		hnennOverDeltaT = hnenn / deltaT;*/
/**/
/*        for ( i = 0; i < n_equations; ++i ){*/
/*			if ( i == 12 ) {*/
/*//				a[k][i][i] += hnennOverDeltaT;*/
/*			}*/
/*//			locDy[k][i] += hnenn * yPrime[k][i];*/
/*		}*/
/*    }*/
/*}*/

/*void TNewton::Update_timedepRHS( TTimePtr time )*/
/*{*/
/*	TGridPtr	locGrid = grid->GetCurrentGrid();	*/
/*	int 		i, k;*/
/*	int			nGridPoints = locGrid->GetNGridPoints();*/
/*	Double 		h;*/
/*	Double 		hm;*/
/*	Double 		hnenn;*/
/*	Double 		*x = locGrid->GetX()->vec;*/
/*	Double		**locDy = dy->mat;*/
/*	Double		**yPrime = time->GetYPrime()->mat;*/
/*	*/
/*    for ( k = 0; k < nGridPoints; ++k ){*/
/*	*/
/*        if ( k == 0 ){*/
/*            hm = x[0] - GetLeft();*/
/*            h = x[1]- x[0];*/
/*        }*/
/*        else if ( k == nGridPoints-1 ){*/
/*            hm = x[k] - x[k-1];*/
/*            h = GetRight() - x[k];*/
/*        }*/
/*        else{*/
/*            hm = x[k] - x[k-1];*/
/*            h = x[k+1] - x[k];*/
/*        }*/
/*		hnenn = h * hm * ( h + hm );*/
/**/
/*        for ( i = 0; i < n_equations; ++i ){*/
/*			locDy[k][i] += hnenn * yPrime[k][i];*/
/*		}*/
/*    }*/
/*}*/

int TNewton::Update_Y( TDampPtr damp, TTimePtr time )
{
    int i, k;
	TGridPtr 	currentGrid = grid->GetCurrentGrid();
	int			currentGridPoints = currentGrid->GetNGridPoints();
	Double		**currentY = currentGrid->GetY()->mat;
	Double		**locDy = dy->mat;
	int			equations = GetNEquations();

	if ( timedep_flag ) {
//		Double	fact;
//		Double	**curryPrime = time->GetYPrime()->mat;
		if ( damp_flag ) {
			Double	lambda = damp->GetLambda();
//			fact = lambda / time->GetDeltaT(); 
			for ( k = 0; k < currentGridPoints; ++k ) {
				for ( i = 0; i < equations; ++i ) {
					currentY[k][i] += lambda * locDy[k][i];
//					curryPrime[k][i] += fact * locDy[k][i];
				}
			}
		}
		else {
//			fact = 1.0 / time->GetDeltaT(); 
			for ( k = 0; k < currentGridPoints; ++k ){
				for ( i = 0; i < equations; ++i ){
					currentY[k][i] += locDy[k][i];
//					curryPrime[k][i] += fact * locDy[k][i];
				}
			}
		}
	}
	else {
		if ( damp_flag ) {
			Double	lambda = damp->GetLambda();
			for ( k = 0; k < currentGridPoints; ++k ) {
				for ( i = 0; i < equations; ++i ) {
					currentY[k][i] += lambda * locDy[k][i];
				}
			}
		}
		else {
			for ( k = 0; k < currentGridPoints; ++k ){
				for ( i = 0; i < equations; ++i ){
					currentY[k][i] += locDy[k][i];
				}
			}
		}
	}

    return 0;
}

void TNewton::BackSolve( MatrixPtr locDy, void *object, Flag updateRHS )
{
        int err = 0;

	if ( updateRHS ) {
		UpdateRHS( locDy, object );
	}

#ifdef NEWSCALEDSOLUTION
	ScaledY( locDy->mat, fJacobiScaler->mat );
#endif

        err = solbt( fA, fB, fC, locDy, ip );
        if ( err ) fprintf( stderr, "# error in solbt\n" );
}

void TNewton::SetNodeInfo( void *object, int k )
{
	TGridPtr	currentGrid = grid->GetCurrentGrid();
	Double		*x = currentGrid->GetX()->vec;
	Double		**y = currentGrid->GetY()->mat;
	Double		**locDy = dy->mat;
	int			N = currentGrid->GetNGridPoints();
	Double		*yLeft = currentGrid->GetYLeft()->vec;
	Double		*yRight = currentGrid->GetYRight()->vec;

	fUtFuncs->SetObjNodeInfo( k, object );

	fNodeInfo->nOfEquations = GetNEquations();
	fNodeInfo->parameter = GetParameter();
	fNodeInfo->gridPoint = k;
	
	if ( k >= 0 && k < N ) {
		fNodeInfo->bcFlagLeft = currentGrid->GetBcFlagLeft();
		fNodeInfo->bcFlagRight = currentGrid->GetBcFlagRight();
		fNodeInfo->bcLeft = currentGrid->GetBcLeft()->vec;
		fNodeInfo->bcRight = currentGrid->GetBcRight()->vec;

		fNodeInfo->a = fA->tensor[k];
		fNodeInfo->b = fB->tensor[k];
		fNodeInfo->c = fC->tensor[k];
		fNodeInfo->rhs = locDy[k];
		fNodeInfo->rhsSaved = dySaved->mat[k];
		fNodeInfo->x = &x[k];
		fNodeInfo->y = y[k];
		fNodeInfo->yMax = fYMax->vec;
	

		if ( k == 0 ) {
			fNodeInfo->firstPoint = TRUE;
			fNodeInfo->lastPoint = FALSE;
			fNodeInfo->h = x[1] - x[0];
			fNodeInfo->hm = x[0] - GetLeft();
			fNodeInfo->hNext = x[2] - x[1];
			fNodeInfo->hmPrev = -1.0e100;
			fNodeInfo->yPrevPrev = NULL;
			fNodeInfo->yPrev = yLeft;
			fNodeInfo->yNext = y[1];
			fNodeInfo->yNextNext = y[2];
			fNodeInfo->bPrev = NULL;
			fNodeInfo->cNext = fC->tensor[1];
		}
		else if ( k == N-1 ) {
			fNodeInfo->firstPoint = FALSE;
			fNodeInfo->lastPoint = TRUE;
			fNodeInfo->h = GetRight() - x[N-1];
			fNodeInfo->hm = x[N-1] - x[N-2];
			fNodeInfo->hNext = -1.0e100;
			fNodeInfo->hmPrev = x[N-2] - x[N-3];
			fNodeInfo->yPrevPrev = y[N-3];
			fNodeInfo->yPrev = y[N-2];
			fNodeInfo->yNext = yRight;
			fNodeInfo->yNextNext = NULL;
			fNodeInfo->bPrev = fB->tensor[N-2];
			fNodeInfo->cNext = NULL;
		}
		else {
			fNodeInfo->firstPoint = FALSE;
			fNodeInfo->lastPoint = FALSE;
			fNodeInfo->h = x[k+1] - x[k];
			fNodeInfo->hm = x[k] - x[k-1];
			if ( k > 1 ) {
				fNodeInfo->hmPrev = x[k-1] - x[k-2];
				fNodeInfo->yPrevPrev = y[k-2];
			}
			if ( k < N-2 ) {
				fNodeInfo->hNext = x[k+2] - x[k+1];
				fNodeInfo->yNextNext = y[k+2];
			}
			fNodeInfo->yPrev = y[k-1];
			fNodeInfo->yNext = y[k+1];
			fNodeInfo->bPrev = fB->tensor[k-1];
			fNodeInfo->cNext = fC->tensor[k+1];
		}
		fNodeInfo->hnenn = fNodeInfo->h * fNodeInfo->hm * ( fNodeInfo->h + fNodeInfo->hm );
	}
}

void TNewton::PrintSolution( Double *x, Double **y, ConstStringArray names, FILE *fp )
{
	TGridPtr	currentGrid = grid->GetCurrentGrid();
    int         i, k;
    int			M = GetNVariables();
    int			N = currentGrid->GetNGridPoints();
	Double		*yleft = currentGrid->GetYLeft()->vec;
	Double		*yright = currentGrid->GetYRight()->vec;
	int 		*len = new int[M];
	char		format[128];

    fprintf( fp, "*\n" );
	
	fprintf( fp, "%-12s", "y" );
	for ( i = 0; i < M; ++i ) {
		len[i] = maxint( strlen( names[i] ), 12 );
		sprintf( format, "\t%%-%ds", len[i] );
//		cerr << "format" << format << NEWL;
		fprintf( fp, format, names[i] );
	}

    fprintf( fp, "\n%-12E", GetLeft() );
    for ( i = 0; i < M; ++i ){
		sprintf( format, "%s%%-%dE", "\t", len[i] );
        fprintf( fp, format, yleft[i] );
    }
    for ( k = 0; k < N; ++k ){
        fprintf( fp, "\n%-12E", x[k] );
        for ( i = 0; i < M; ++i ){
			sprintf( format, "%s%%-%dE", "\t", len[i] );
            fprintf( fp, format, y[k][i] );
        }
    }
    fprintf( fp, "\n%-12E", GetRight() );
    for (i = 0; i < M; ++i ){
		sprintf( format, "%s%%-%dE", "\t", len[i] );
        fprintf( fp, format, yright[i] );
    }
    fprintf( fp, "\n" );
	
	delete len;
}

void TNewton::PrintResiduals( MatrixPtr rhs, ConstStringArray names )
{
    int			i, k, rows = rhs->rows, cols = rhs->cols;
    Double		norm = 0.0;
	Double		**mat = rhs->mat;
	FILE		*fp;
	static int	init = FALSE;

	if ( !init ) {
		sprintf( GetOutFileBuff(), "%sresiduals.tout", GetOutputPath() );
		fp = fopen( GetOutFileBuff(), "w" );
		
		fprintf( fp, "%-4s", "iter" );
		for ( i = 0; i < rows; ++i ) {
			fprintf( fp, "\t%-12s", names[i] );
		}
		fprintf( fp, "\n" );
		init = TRUE;
	}
	else {
		sprintf( GetOutFileBuff(), "%sresiduals.tout", GetOutputPath() );
		fp = fopen( GetOutFileBuff(), "a" );
	}

	fprintf( fp, "%-4d", n_iter );
	for ( i = 0; i < rows; ++i ) {
		for ( k = 0, norm = 0.0; k < cols; ++k ) {
            norm += MIN(fabs(mat[k][i]), 1.0e99) * MIN(fabs(mat[k][i]), 1.0e99);
//            norm += mat[k][i] * mat[k][i];
		}
		norm /= cols;
		norm = sqrt( norm );
		fprintf( fp, "\t%-12E", norm );
    }
	fprintf( fp, "\n" );
	
	fclose( fp );
}

void TNewton::ScaleMatrix( MatrixPtr matPtr, VectorPtr vecPtr )
{
	int		i, k;
	int		gridPoints = matPtr->cols;
	int		variables = matPtr->rows;
	Double	**y = matPtr->mat;
	Double	*scaler = vecPtr->vec;

	if ( vecPtr->len < variables ) {
		NewtonFatalError( "invalid scaler for current matrix" );	
	}
	
	for ( i = 0; i < variables; ++i ) {
		for ( k = 0; k < gridPoints; ++k ) {
			y[k][i] /= scaler[i];
		}
	}
}

void TNewton::DeScaleMatrix( MatrixPtr matPtr, VectorPtr vecPtr )
{
	int		i, k;
	int		gridPoints = matPtr->cols;
	int		variables = matPtr->rows;
	Double	**y = matPtr->mat;
	Double	*scaler = vecPtr->vec;

	if ( vecPtr->len < variables ) {
		NewtonFatalError( "invalid scaler for current matrix" );	
	}
	
	for ( i = 0; i < variables; ++i ) {
		for ( k = 0; k < gridPoints; ++k ) {
			y[k][i] *= scaler[i];
		}
	}
}

void TNewton::SetNOfEquations( int number )
{
	if ( number > n_equations ) {
		NewtonFatalError( "invalid call of function SetNOfEquations" );	
	}
	else {
		fA->rows = number;
		fA->cols = number;
		fB->rows = number;
		fB->cols = number;
		fC->rows = number;
		fC->cols = number;
		dy->rows = number;
	}
}

void TNewton::AdjustNewtonDimensions( int nOfPoints )
{
	dy->cols = nOfPoints;
	fA->planes = nOfPoints;
	fB->planes = nOfPoints;
	fC->planes = nOfPoints;
}

void TNewton::ScaleRHS( MatrixPtr dys )
{
	int		i, k;
	int		nEquations = n_equations;
	int		nGridPoints = GetCurrentGridPoints();
	Double	**rhs = dys->mat;
	Double	**scaler = fJacobiScaler->mat;

	for ( k = 0; k < nGridPoints; ++k ) {
		for ( i = 0; i < nEquations; ++i ) {
			rhs[k][i] *= scaler[k][i];
		}
	}
}

void TNewton::NewDeScaleSolution( Double **dy )
{
	int		i, k;
	int		nEquations = fA->rows;
	int		nGridPoints = fA->planes;
	Double	**scaler = fJacobiScaler->mat;
	
	for ( k = 0; k < nGridPoints; ++k ) {
		for ( i = 0; i < nEquations; ++i ) {
			dy[k][i] *= scaler[k][i];
		}
	}
}

void TNewton::ScaleSystem( void )
{
	int		i, j, k;
	int		nEquations = fA->rows;
	int		nGridPoints = fA->planes;
	Double	***a = fA->tensor;	//	a[gridPoint][variable][equation]
	Double	***b = fB->tensor;
	Double	***c = fC->tensor;
	Double	**dY = dy->mat;
	Double	**scaler = fJacobiScaler->mat;
	
	for ( k = 0; k < nGridPoints; ++k ) {
		for ( i = 0; i < nEquations; ++i ) {	// loop equation
			scaler[k][i] = ( fabs(a[k][i][i]) > 1e-15 ) ? 1.0 / a[k][i][i] : 1.0;
			dY[k][i] *= scaler[k][i];
			for ( j = 0; j < nEquations; ++j ) {	// loop variable
				a[k][j][i] *= scaler[k][i];
				b[k][j][i] *= scaler[k][i];
				c[k][j][i] *= scaler[k][i];
			}
		}
	}
}

void TNewton::ScaledY( Double **dY, Double **scaler )
{
	int		i, k;
	int		nEquations = fA->rows;
	int		nGridPoints = fA->planes;
	
	for ( k = 0; k < nGridPoints; ++k ) {
		for ( i = 0; i < nEquations; ++i ) {	// loop equation
			dY[k][i] *= scaler[k][i];
		}
	}
}

void TNewton::TestSolution( void )
{
	int		i, j, k;
	int		nEquations = fA->rows;
	int		nGridPoints = fA->planes;
	int 	lastPoint = nGridPoints - 1;
	Double	***a = fA->tensor;	//	a[gridPoint][variable][equation]
	Double	***b = fB->tensor;
	Double	***c = fC->tensor;
	//	Double	**rhs = dy->mat;
	Double	**mat = grid->GetWorkMatrix()->mat;
	Double	**y = grid->GetCurrentGrid()->GetY()->mat;
	
	ClearMatrix( grid->GetWorkMatrix() );
	
	for ( k = 0; k < nGridPoints; ++k ) {
		for ( i = 0; i < nEquations; ++i ) {	// loop equation
			for ( j = 0; j < nEquations; ++j ) {	// loop variable
				mat[k][i] += y[k][j] * a[k][j][i];
				if ( k != lastPoint ) {
					mat[k][i] += y[k+1][j] * b[k][j][i];
				}
				if ( k != 0 ) {
					mat[k][i] += y[k-1][j] * c[k][j][i];
				}
			}
		}
	}
}

void TNewton::PrintTheVector( VectorPtr vecPtr, char *header )
{
	Double		*vec = vecPtr->vec;
	TGridPtr	currentGrid = grid->GetCurrentGrid();
	Double		*x = currentGrid->GetX()->vec;
	int			gridPoints = currentGrid->GetNGridPoints();
	FILE		*fp = NULL;
	
	if ( gridPoints > vecPtr->len ) {
		fprintf( stderr, "#warning: vector not printed, because it's too short for the current grid\n" );
		return;
	}
	if ( strlen( header ) > 27 ) {
		NewtonFatalError( "invalid filename" );
	}
	sprintf( GetOutFileBuff(), "%s%s.dout", GetOutputPath(), header );
	fp = fopen( GetOutFileBuff(), "w" );

	fprintf( fp, "*\n" );
	fprintf( fp, "eta\t%s\n", header );
	for ( int k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "%-.6e\t%-.6e\n", x[k], vec[k] );
	}
	
	fclose( fp );
}

void TNewton::PrintVectorOfMatrix( MatrixPtr matPtr, int number, char *header )
{
	Double		**mat = matPtr->mat;
	TGridPtr	currentGrid = grid->GetCurrentGrid();
	Double		*x = currentGrid->GetX()->vec;
	int			gridPoints = currentGrid->GetNGridPoints();
	FILE		*fp = NULL;
	
	if ( gridPoints > matPtr->cols ) {
		fprintf( stderr, "#warning: vector not printed, because it's too short for the current grid\n" );
		return;
	}
	if ( strlen( header ) > 27 ) {
		NewtonFatalError( "invalid filename" );
	}
	sprintf( GetOutFileBuff(), "%s%s.dout", GetOutputPath(), header );
	fp = fopen( GetOutFileBuff(), "w" );

	fprintf( fp, "*\n" );
	fprintf( fp, "eta\t%s\n", header );
	for ( int k = 0; k < gridPoints; ++k ) {
		fprintf( fp, "%-.6e\t%-.6e\n", x[k], mat[k][number] );
	}
	
	fclose( fp );
}

void TNewton::CheckGridding( TTimePtr tim )
{
	if ( GetNIter() == 1 || GetNIter() == GetMaxIter() ) {
		fIsGriddingStep = FALSE;
//		grid->UnSetSolHasNewGrid();
	}
	else if ( ( GetNIter() - grid->GetIterOfGeneration() - GetDeltaNewGrid() ) >= 0 ) {
		fIsGriddingStep = TRUE;
	}
	else if ( grid->GetOneSolOneGrid() && !grid->GetSolHasNewGrid() 
				&& ( ( normRes < 1.0e-1 && normDy < 1.0e-1 ) || normDy < 100.0 * TolDy ) 
				&& ( ( GetTimedepFlag() ) ? ( tim->GetTimeConverged() 
				|| tim->GetTimeStep() % 10  == 0 
				|| tim->GetDeltaT() >= tim->GetDeltaTMax() ) : TRUE )
				) {
	//		cerr << "third:" << grid->GetOneSolOneGrid() << TAB << grid->GetSolHasNewGrid() << TAB << normRes << TAB << normDy << NEWL;
		fIsGriddingStep = TRUE;
	}
	else {
		fIsGriddingStep = FALSE;
	}
}

void TNewton::WriteOutput( void *object, FILE *fp, char* tail )
{
	fUtFuncs->WriteOutput( object, fp, tail );
}

void TNewton::SetUtFuncs( JacFuncPtr jFirst, JacFuncPtr jRest, JacFuncPtr jLast
					, RHSFuncPtr rFirst, RHSFuncPtr rRest, RHSFuncPtr rLast
					, OutputFuncPtr writeOutput, PostIterFuncPtr postIter
					, SetObjNodeInfoPtr setObjNodeInfo, PostConvFuncPtr postConvergence
					, GetVariableNamesPtr getVariableNames
					, UpdateBoundaryPtr updateLeftBoundary
					, UpdateBoundaryPtr updateRightBoundary )
{
	fUtFuncs->JacobiFirstPoint = jFirst;
	fUtFuncs->JacobiRestPoint = jRest;
	fUtFuncs->JacobiLastPoint = jLast;

	fUtFuncs->RHSFirstPoint = rFirst;
	fUtFuncs->RHSRestPoint = rRest;
	fUtFuncs->RHSLastPoint = rLast;
	
	fUtFuncs->WriteOutput = writeOutput;
	fUtFuncs->PostIter = postIter;
	fUtFuncs->SetObjNodeInfo = setObjNodeInfo;
	fUtFuncs->PostConvergence = postConvergence;
	fUtFuncs->UpdateLeftBoundary = updateLeftBoundary;
	fUtFuncs->UpdateRightBoundary = updateRightBoundary;
	fUtFuncs->GetVariableNames = getVariableNames;
}

ConstStringArray TNewton::GetVariableNames( void *object )
{
	return fUtFuncs->GetVariableNames( object );
}

FILE *TNewton::GetOutfile( const char *name, FileType type, Flag iter )
{
	char	*fullName = GetOutfileName( name, type, iter );
	FILE	*fp;

	if ( !( fp = fopen( fullName, "w") ) ) { 
		cerr << "#warning: unable to open file " << fullName << NEWL;
		exit(2);
	}
	
	return fp;
}

char *TNewton::GetOutfileName( const char *name, FileType type, Flag iter )
{
#if defined (applec) || defined (powerc)
	Double	dummy;
#else
	Double		dummy;
#endif
	int			counter;
	
	if ( iter ) {
		counter = ( int ) ( modf( GetNIter()/10.0, &dummy ) * 10.0 );
		sprintf( GetOutFileBuff(), "%s%s%d.%cout"
				, GetOutputPath(), name, counter, ( type == kData ) ? 'd' : 't' );
	}
	else {
		sprintf( GetOutFileBuff(), "%s%s.%cout"
				, GetOutputPath(), name, ( type == kData ) ? 'd' : 't' );
	}
	
	return GetOutFileBuff();
}

void TBVPSolver::ReInit( void )
{
	bt->UnSetLeaveContin();
	bt->InitNIter();
	bt->GetGrid()->UnSetSolHasNewGrid();
	if ( bt->GetTimedepFlag() ) {
		time->UnSetTimeConverged();
		time->InitDeltaT();
		time->InitNTimeStep();

	}
}

void TNewton::SetMaximaOfSolution( void )
{
	MatrixPtr		YMat = grid->GetCurrentGrid()->GetY();
	Double			**Y = YMat->mat;
	Double			*locY;
	Double			*YMax = fYMax->vec;
	
	locY = Y[0];
	for ( int j = 0; j < YMat->rows; ++j ) {
		YMax[j] = locY[j];
	}
	for ( int k = 1; k < YMat->cols; ++k ) {
		locY = Y[k];
		for ( int j = 0; j < dy->rows; ++j ) {
			if ( locY[j] > YMax[j] ) YMax[j] = locY[j];
		}
	}
}
