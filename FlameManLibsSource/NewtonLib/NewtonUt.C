#include "Newton.h"

#define NEWMONITOR

#define NEWMAXCHANGE

#define PREDICTOR

UTFuncs::UTFuncs( void )
{
	JacobiFirstPoint = NULL;
	JacobiRestPoint = NULL;
	JacobiLastPoint = NULL;

	RHSFirstPoint = NULL;
	RHSRestPoint = NULL;
	RHSLastPoint = NULL;

	WriteOutput = NULL;
	PostIter = NULL;
	SetObjNodeInfo = NULL;
	PostConvergence = NULL;
	GetVariableNames = NULL;
}

void TDamp::InitDamp( TNewtonPtr bt, TBVPInputPtr input )
{
	int	maxGridPoints = bt->GetMaxGridPoints();
	int variables = bt->GetNVariables();

	Lambda = 1;
    LambdaMin = input->fLambdaMin;
    Sigma = 0.01;
    Tau = 0.1;
	converg_Lambda = TRUE;
    yOld = NewMatrix( variables, maxGridPoints, kColumnPointers );
    dyOld = NewMatrix( variables, maxGridPoints, kColumnPointers );
}

TDamp::~TDamp( void )
{
	DisposeMatrix( dyOld );
	DisposeMatrix( yOld );
}

int TDamp::Update_Lambda( TNewtonPtr bt, TTimePtr time, void *object )
{
//    int                 i, k;
        Double          Lambda_new, crit_func0, crit_func, crit_funcMax;
	TGridPtr	fineGrid = bt->GetGrid()->GetCurrentGrid();
	int             M = bt->GetNEquations();
	int             N = fineGrid->GetNGridPoints();
	VectorPtr	scaler = bt->GetGrid()->GetSolutionScaler();
	Double		**yFine = fineGrid->GetY()->mat;
	Double		**dy = bt->GetDy()->mat;
	
	copy_mat( yOld, fineGrid->GetY() );
	if ( bt->GetNIter() == 1 || bt->isNewGrid() ){
	        Lambda = 0.1;
/*		if ( bt->GetTimedepFlag() ) {*/
/*			fprintf( stderr, "attention: damping and time dependence should not be used together\n" );*/
/*		}*/
		bt->UpdateSolution( object, this, time );
		bt->BackSolve( dyOld, object );
	}
	else if ( !bt->fullstep ) {
	        Lambda = 1.0;
		bt->UpdateSolution( object, this, time );
		bt->BackSolve( dyOld, object );
	}
	else {

		crit_func0 = ScaledNorm( bt->GetDy(), scaler );
		crit_func0 *= 0.5 * crit_func0;

		Predict_Lambda( bt, bt->GetDy(), bt->GetDyNorm() );
		//cerr << "predicted Lambda is " << Lambda << NEWL; //GB
		for (;;){
			bt->UpdateSolution( object, this, time );
			
			//  compute crit_func( y + Lambda*dy )
			bt->BackSolve( dyOld, object );

			crit_func = ScaledNorm( dyOld, scaler );
			crit_func *= crit_func * 0.5;

			crit_funcMax = ( 1.0 - 2.0 * Lambda * Sigma ) * crit_func0;
			if ( crit_func > crit_funcMax ){
			       Lambda_new = Lambda * Lambda * crit_func0 /
				 ( ( 2.0 * Lambda - 1.0 ) * crit_func0 + crit_func );
			       Lambda = ( Tau * Lambda > Lambda_new ) ? Tau * Lambda : Lambda_new;
			       if ( Lambda < LambdaMin ){
				     //converg_Lambda = 0;
				     //                   fprintf ( stderr, "# no convergence\n" );
				     Lambda = 10.0 * LambdaMin;
				     break;
			       }
			}
			else{
			       break;  //  Lambda is good                                      
			}
			}
	}
    copy_mat( fineGrid->GetY(), yOld );
    return 0;
}

void TDamp::Predict_Lambda( TNewtonPtr bt, MatrixPtr dy, Double norm )
{
    Double mu;

    sub_mat( 2, dy, dyOld );

    mu = norm / ScaledNorm( dyOld, bt->GetGrid()->GetSolutionScaler() );
    Lambda = MAX( Tau, mu ) * Lambda;
    Lambda = MIN( Lambda, 1.0 );
    Lambda = MAX( Lambda, LambdaMin );
}

void TDamp::AdjustDampDimension( int nOfPoints )
{
	yOld->cols = nOfPoints;
	dyOld->cols = nOfPoints;
}

void TTime::InitTime( TNewtonPtr bt, TBVPInputPtr input )
{
/*  don't change the following parameters  */
	int	maxGridPoints = bt->GetMaxGridPoints();
	int variables = bt->GetNVariables();

    stepsGood = 0;
    stepsBad = 0;
    Delta_t = fDeltaTStart = input->fDeltaTStart;
    normGoUp = 1000.0;
    normOld = 1000.0;
    Delta_tMin = MIN( 1.0e-10, Delta_t );
    Delta_tMax = 10.0;
	fNTimeStep = 1;
	fAllowedFact = 2.0;
	InitDeltaTMax( bt );
	fIsDtGuess = FALSE;
	fTimeConverged = FALSE;
	UnsetAdjustTime();    

	yOld = NewMatrix( variables, maxGridPoints, kColumnPointers );
    yPrime = NewMatrix( variables, maxGridPoints, kColumnPointers );
    x_old = NewVector( maxGridPoints );
}

TTime::~TTime( void )
{
	DisposeVector( x_old );
	DisposeMatrix( yOld );
	DisposeMatrix( yPrime );
}

/*void TTime::InitYPrime( TNewtonPtr bt, void *object )*/
/*{*/
/*//    bt->UpdateRHS( yPrime, object );*/
/*	ClearMatrix( yPrime );*/
/*	MatrixPtr tmpMat = NewMatrix( bt->GetNVariables(), bt->GetMaxGridPoints(), kColumnPointers );*/
/*    bt->UpdateRHS( tmpMat, object );*/
/*	copy_mat( yPrime, tmpMat );*/
/*	for ( int k = 0; k < bt->GetCurrentGridPoints(); ++k ){*/
/*		bt->SetNodeInfo( object, k );*/
/*		for ( int i = 0; i < bt->GetNVariables(); ++i ) {*/
/*			yPrime->mat[k][i] /= -bt->GetNodeInfo()->hnenn;*/
/*		}*/
/*	}*/
/*	int 	co, ro;*/
/*	Double	locNorm = MaxNorm_mat( yPrime, &co, &ro );*/
/*	*/
/*	fprintf( stderr, "yPrime: norm = %g at yPrime[%d][%d]\n"*/
/*				, locNorm, co, ro );*/
/*}*/

void TTime::AdjustTimeStep( TBVPSolverPtr solver, Flag converged, void *object )
{
	TNewtonPtr	bt = solver->bt;
	TGridPtr	fineGrid = bt->GetGrid()->GetCurrentGrid();
//	Double		locNorm;
	Double		dtguess;
	Double		maxChange;
    if ( !converged ) {
		fprintf( stderr, "to old solution\n" );
		bt->GetGrid()->ToOldSolution( solver, x_old, yOld );
		bt->GetUTFuncs()->PostIter( object );
		bt->WriteOutput( object, NULL, "_old" );
		bt->GetGrid()->UnSetSolHasNewGrid();

/*		FILE	*fpr = NULL;*/
/*		sprintf( bt->GetOutFileBuff(), "%sNewSolution.dout", bt->GetOutputPath() );*/
/*		if ( ( fpr = fopen( bt->GetOutFileBuff(), "w" ) ) == NULL ) {*/
/*			cerr << "#warning: can't open file NewSolution.dout" << NEWL;*/
/*		}*/
/*		bt->PrintSolution( fineGrid->GetX()->vec, fineGrid->GetY()->mat, bt->GetVariableNames( object ), fpr );*/
/*		fclose( fpr );*/

		dtguess = MIN( 0.1 * Delta_t, 1.0 );
		if ( Delta_t < Delta_tMin ) {
		  fprintf( stderr, "Delta_t = %g < Delta_tMin = %g\n", Delta_t, Delta_tMin );
		  exit(2);
		}
		fprintf( stderr, "\nrepeat time step no. %d with dt = %g\n", fNTimeStep, dtguess );
	
		Delta_tOldOld = Delta_tOld;
		Delta_tOld = Delta_t;

		Delta_t = dtguess;
		if ( fIsDtGuess ) {
			fAllowedFact = MAX( 1.1, 0.9 * fAllowedFact );
		}
		fIsDtGuess = FALSE;

		ReInitTimeMats( bt, object, Delta_t, Delta_tOld );
	}
    else {
		MatrixPtr	locOlddeb = NewMatrix( yOld->phys_rows, yOld->phys_cols, kColumnPointers );
		copy_mat( locOlddeb, yOld );
		fprintf( stderr, "timestep converged\n" );
		bt->WriteOutput( object, NULL, "_t" );
		int co, ro;
/*		locNorm = MaxNorm_mat( yPrime, &co, &ro );*/
/*				char		name[128];*/
/*				sprintf( name, "YPrime.tout" );*/
/*	*/
/*				AMPrintOptions	prnt;*/
/*				DefaultAMPOpts( &prnt );*/
/*				FILE *fp = fopen( name, "w" );*/
/*				fprintf( fp, "print YPrime\n" );*/
/*				PrintMatrix( yPrime, &prnt, fp );*/
/*				fclose( fp );*/
/*	*/
/*		dtguess = fabs( fineGrid->GetY()->mat[co][ro] / yPrime->mat[co][ro] );*/
/*		fprintf( stderr, "yPrime: norm = %g at yPrime[%d][%s]\tdt_guess = %g\n"*/
/*				, locNorm, co, bt->GetVariableNames( object )[ro], dtguess );*/
		
		// save solution
/*				sprintf( name, "YOverYOld.tout" );*/
/*	*/
/*				fp = fopen( name, "w" );*/
/*				fprintf( fp, "print Y\n" );*/
/*				PrintMatrix( fineGrid->GetY(), &prnt, fp );*/
/*				fprintf( fp, "print Yold\n" );*/
/*				PrintMatrix( yOld, &prnt, fp );*/
//		div_mat( 2 ,fineGrid->GetY(), yOld );
		maxChange = MaxChange( yOld, fineGrid->GetY(), &co, &ro );
//		dtguess = Delta_t * 2.0 / maxChange;
//		dtguess = Delta_t / ( maxChange - 1.0 );
//fprintf( stderr, "dtguess2 = %g\n", dtguess );

/*				fprintf( fp, "print Y/Yold\n" );*/
/*				PrintMatrix( yOld, &prnt, fp );*/
/*				fclose( fp );*/
//		dtguess = maxChange / ( maxChange - 1.0 ) * Delta_t;
		dtguess = ( fAllowedFact - 1.0 ) * Delta_t / MAX( 1.0e-10, maxChange - 1.0 );
		fprintf( stderr, "Y/YOld = %g at Y[%d][%s]  Y = %g  YOld = %g  dtguess = %g\n"
						, maxChange, co, bt->GetVariableNames( object )[ro]
						, fineGrid->GetY()->mat[co][ro]
						, locOlddeb->mat[co][ro], dtguess );

// set new time step
/*		if ( bt->GetNIter() > 10 ) {*/
/*			Delta_t *= 0.5;*/
/*		}*/
/*		else {*/
/*			if ( maxChange < 1.01 ) {*/
/*				Delta_t *= 4.0;*/
/*			}*/
/*			else if ( bt->GetNIter() < 4 ) {*/
/*				Delta_t *= 2.0;*/
/*			}*/
/*			else if ( maxChange < 1.4 ) {*/
/*				Delta_t *= 2.0;*/
/*			}*/
/*			else if ( maxChange > 4.0 ) {*/
/*				Delta_t *= 0.5;*/
/*			}*/
/*//			Delta_t = dtguess;*/
/*		}*/
//		Delta_t = 0.01 * dtguess;
		DisposeMatrix( locOlddeb );
		++fNTimeStep;

//		if ( maxChange < 1.1 ) {
//			Delta_t = MIN( 5.0 * Delta_t, dtguess);
//		}
//      else
		if ( maxChange < 2.0 ) {
		  if (maxChange < 1.01 ) {
		    Delta_t = MIN( 5.0 * Delta_t, dtguess);
		  } 
		  else {
		    Delta_t = MIN( 2.0 * Delta_t, dtguess);
		  }
		}
		// else Delta_t = Delta_t;

		//		if ( Delta_t >= Delta_tMax && !GetTimeConverged() && solver->bt->GetGrid()->GetSolHasNewGrid() ) {

		if ( !GetTimeConverged() ) {
		  if ( Delta_t >= Delta_tMax && (bt->GetGrid()->GetOneSolOneGrid()?(bt->GetGrid()->GetSolHasNewGrid()):TRUE) ) {
		    //if ( solver->bt->GetGrid()->GetSolHasNewGrid() ) {
//		    if ( bt->GetGrid()->GetOneSolOneGrid()?(bt->GetGrid()->GetSolHasNewGrid()):TRUE ) {
		      cerr << NEWL << "switch off time dependence" << NEWL;
		      bt->GetGrid()->UnSetSolHasNewGrid();
		      SetTimeConverged();
		      ReInitTimeMats( bt, object, Delta_t, Delta_tOld );
		      if ( solver->damp ) bt->SetDampFlag( TRUE ) ;
//		    }
		  }
		  else {
		  fprintf(stderr, "step = %d\nmod = %d\n", GetTimeStep(), GetTimeStep() % 10);
  		  	if ( ( GetTimeStep()) % 10 == 0 ) { 
				bt->GetUTFuncs()->JacobiFirstPoint( object, NULL );	
			}
		    fprintf( stderr, "\nstart time step no. %d\n", fNTimeStep );
		    ReInitTimeMats( bt, object, Delta_t, Delta_tOld );
		  }
		}

		fIsDtGuess = TRUE;
	
		Delta_tOldOld = Delta_tOld;
		Delta_tOld = Delta_t;
	}


//	if ( !GetTimeConverged() ) {
//		ReInitTimeMats( bt, object, Delta_t, Delta_tOld );
//	}
}

void TTime::ReInitTimeMats( TNewtonPtr bt, void *object, Double dt, Double dtOld )
{
	bt->InitNIter();
	//if ( fmod( fNTimeStep, 5 ) == 0 ) bt->GetGrid()->UnSetSolHasNewGrid();
	if ( ( fNTimeStep % 5 ) == 0 ) bt->GetGrid()->UnSetSolHasNewGrid();
	
	copy_mat( yPrime, yOld );
	copy_mat( yOld, bt->GetGrid()->GetCurrentGrid()->GetY() );
	copy_vec( x_old, bt->GetGrid()->GetCurrentGrid()->GetX() );

//			fprintf( stderr, "save solution\n" );
			FILE	*fpr = NULL;
			sprintf( bt->GetOutFileBuff(), "%sSavedsolution.dout", bt->GetOutputPath() );
			if ( ( fpr = fopen( bt->GetOutFileBuff(), "w" ) ) == NULL ) {
				cerr << "#warning: can't open file Savedsolution.dout" << NEWL;
			}
			Double	*x = x_old->vec;
			bt->PrintSolution( x, yOld->mat, bt->GetVariableNames( object ), fpr );
			fclose( fpr );

#ifdef PREDICTOR
    int			i, k, rows = yPrime->rows, cols = yPrime->cols;
	TGridPtr	fineGrid = bt->GetGrid()->GetCurrentGrid();
	Double		**y = fineGrid->GetY()->mat;
	Double		**y0 = yPrime->mat;
	Double		delta_tOverDeltatOld = dt / dtOld;

	if ( fNTimeStep != 2 ) { // fNTimeStep has already been incremented
	for ( i = 0; i < rows; ++i ) {
		for ( k = 0; k < cols; ++k ) {
			y[k][i] = y[k][i] + delta_tOverDeltatOld * ( y[k][i] - y0[k][i] ); /* yPrime contains yOld */
		}
	}
	}
	
	bt->GetUTFuncs()->PostIter( object );
#else
	ClearMatrix( yPrime );
#endif
}

void TTime::InitTimeMats( TNewtonPtr bt, void */*object*/ )
{
	copy_mat( yOld, bt->GetGrid()->GetCurrentGrid()->GetY() );
	copy_vec( x_old, bt->GetGrid()->GetCurrentGrid()->GetX() );

	ClearMatrix( yPrime );
}

/*Double TTime::UpdateYPrime( TNewtonPtr bt )*/
/*{*/
/*	fprintf( stderr, "ERROR: function UpdateYPrime shouldn't be called\n" );*/
/*	TGridPtr 	currentGrid = bt->GetGrid()->GetCurrentGrid();*/
/*	MatrixPtr	locY = currentGrid->GetY();*/
/*	Double		norm;*/
/*	static		Flag init = FALSE;*/
/**/
/*			char		name[128];*/
/*			sprintf( name, "YPrime.tout" );*/
/**/
/*			AMPrintOptions	prnt;*/
/*			DefaultAMPOpts( &prnt );*/
/*			FILE *fp = fopen( name, "w" );*/
/*	*/
/*	if ( !init ) {*/
/*		copy_mat( yOld, locY );*/
/*		init = TRUE;*/
/*	}*/
/*	*/
/*	copy_mat( yPrime, locY );*/
/*			fprintf( fp, "print Y\n" );*/
/*			PrintMatrix( yPrime, &prnt, fp );*/
/*	sub_mat( 1 ,yPrime, yOld );*/
/*			fprintf( fp, "print Yold\n" );*/
/*			PrintMatrix( yOld, &prnt, fp );*/
/*			fprintf( fp, "print Y-Yold\n" );*/
/*			PrintMatrix( yPrime, &prnt, fp );*/
/*			fprintf( fp, "Delta_t = %g\n", Delta_t );*/
/*	mult_mat( yPrime, 1.0 / Delta_t );*/
/*	*/
/*	int co, ro;*/
/*	norm = MaxNorm_mat( yPrime, &co, &ro );*/
/*			fprintf( fp, "print YPrime\n" );*/
/*			PrintMatrix( yPrime, &prnt, fp );*/
/*			fclose( fp );*/
/*	*/
/*	fprintf( stderr, "set yPrime: norm = %g at yPrime[%d][%d]\tdt_guess = %g\n"*/
/*				, norm, co, ro, locY->mat[co][ro] / yPrime->mat[co][ro]  );*/
/*	return norm;*/
/*}*/

void TTime::AdjustTimeDimension( int nOfPoints )
{
	yPrime->cols = nOfPoints;
	yOld->cols = nOfPoints;
	x_old->len = nOfPoints;
}

void TContinuation::InitContin( TNewtonPtr bt, TBVPInputPtr input )
{
	int	maxGridPoints = bt->GetMaxGridPoints();
	int variables = bt->GetNVariables();

    cont_steps = input->fContSteps;

    n_cont = 0;
    inc = 1.0 / cont_steps;
	parameter = 1.0;
	parameter = inc;

    x_last_cont = NewVector( maxGridPoints );
    y_last_cont = NewMatrix( variables, maxGridPoints, kColumnPointers );
	bt->SetParameterPtr( &parameter );
}

TContinuation::~TContinuation( void )
{
	DisposeMatrix( y_last_cont );
	DisposeVector( x_last_cont );
}

int TContinuation::Contin_step( TBVPSolverPtr solver, TDampPtr damp )
{
	TNewtonPtr	bt = solver->bt;
	TGridPtr	fineGrid = bt->GetGrid()->GetCurrentGrid();
    
    if ( n_cont == 0 ){
        ++n_cont;
    }
    else{
        if ( bt->GetConvergNewton() == 1 ){
            if ( bt->GetNIter() < 10 ){
                inc *= 1.1;
            }
			parameter += inc;
			parameter = MIN( 1.0, parameter );
            ++n_cont;
			copy_vec( x_last_cont, fineGrid->GetX() );
            copy_mat( y_last_cont, fineGrid->GetY() );
			solver->ReInit();
        }
        else if ( n_cont == 1){
            /*  if there is no convergence in the first step, set
                the parameter to zero and solve linear system               */
			parameter = 0.0;
        }
        else{
            if ( bt->GetDampFlag() ){
                damp->SetConvergence();
            }
            
			TAdaptiveGridPtr	adapGrid = bt->GetGrid();
			
			adapGrid->ToOldSolution( solver, x_last_cont, y_last_cont );
            inc *= 0.5;
            parameter -= inc;
        	++n_cont;
        }
    }

    fprintf( stderr, "n_cont = %d\n", n_cont );
    bt->InitNIter();

    return 0;
}

void TContinuation::AdjustContinDimension( int nOfPoints )
{
	x_last_cont->len = nOfPoints;
	y_last_cont->cols = nOfPoints;
}

void TBVPInput::InitTBVPInput( int nVariables )
{
	fWriteBT = FALSE;
	fWriteResiduum = FALSE;
	fWatchGridding = FALSE;
	fWriteEverySolution = FALSE;
	fOutputPath = FALSE;

	fWriteFullRes = FALSE;
	fUseModifiedNewton = FALSE;

	bcFlagLeft = new int[ nVariables ];
	bcFlagRight = new int[ nVariables ];
	fBcLeft = NewVector( nVariables );
	fBcRight = NewVector( nVariables );
	yleft = NewVector( nVariables );
	yright = NewVector( nVariables );

	for ( int i = 0; i < nVariables; ++i ) {
		bcFlagLeft[i] = 0;
		bcFlagRight[i] = 0;
		fBcLeft->vec[i] = 0.0;
		fBcRight->vec[i] = 0.0;
		yleft->vec[i] = 0.0;
		yright->vec[i] = 0.0;
	}

	// TNewton
	fNVariables = -1;
	fInitialEquations = -1;
	fMaxGridPoints = 301;
	fInitialGridPoints = 101;
	fDampFlag = FALSE;
	fTimeDepFlag = FALSE;
	fContinFlag = FALSE;
	fDeltaNewGrid = 3;
	fTolRes = 1.0e-15;
    fTolDy = 1.0e-10;
    fMaxIter = 100;
    fLeft = 0.0;
    fRight = 1.0;

	// TAdaptiveGrid
	fOneSolOneGrid = FALSE;
    fR = 60.0;
    fQ = 0.35;

	// TDamp
    fLambdaMin = 0.01;

	// TContinuation
    fContSteps = 10;
}

TBVPInput::~TBVPInput( void )
{
	DisposeVector( yright );
	DisposeVector( yleft );
	DisposeVector( fBcRight );
	DisposeVector( fBcLeft );
	
	delete bcFlagRight;
	delete bcFlagLeft;
}

char *GetFullPath( const char *relPath, Flag type ) 
{
	char		*fullPath;
	char		*home = NULL;
#	if defined (applec) || defined (powerc)
#		ifdef DEBUGPOWERC
		fprintf( stderr, "powerc defined\n" );
		fprintf( stderr, "relpath = %s\n", relPath );
		fprintf( stderr, "type = %d\n", type );
#		endif
		const char	separator = ':';
		const char	wrongSep = '/';
		char work[2];
#	else
		const char	separator = '/';
		const char	wrongSep = ':';
#	endif
//	error checking
	if ( !relPath ) {
		return NULL;
	}
	int len = strlen( relPath );
	int	homeLen = 0;

	if ( strncmp( relPath, "$HOME/", 6 ) == 0 
			|| strncmp( relPath, "$HOME:", 6 ) == 0 ) {
		relPath += 6;
#		ifdef DEBUGPOWERC
		fprintf( stderr, "$HOME found\n" );
#		endif
#	if defined (applec) || defined (powerc)

#	else
		if ( ( home = getenv( "HOME" ) ) == NULL ) NewtonFatalError( "can't find environment variable $HOME" );
		homeLen = strlen( home ) + 1; // +1 for seperator
#	endif
	}
	else if ( strncmp( relPath, "~/", 2 ) == 0 
			|| strncmp( relPath, "~:", 2 ) == 0 ) {
#		ifdef DEBUGPOWERC
		fprintf( stderr, "~: found\n" );
#		endif
		relPath += 2;
#	if defined (applec) || defined (powerc)

#	else
		if ( ( home = getenv( "HOME" ) ) == NULL ) NewtonFatalError( "can't find environment variable $HOME" );
		homeLen = strlen( home ) + 1; // +1 for seperator
#	endif
	}
#	if defined (applec) || defined (powerc)
	else if (  strncmp( relPath, "./", 2 ) == 0 ) {
#		ifdef DEBUGPOWERC
		fprintf( stderr, "./ found\n" );
#		endif
		relPath += 2;
		work[0] = ':';
		work[1] = '\0';
		home = work;
	}
#	endif

	if ( relPath[ strlen( relPath ) - 1 ] == separator 
		|| relPath[ strlen( relPath ) - 1 ] == wrongSep 
		|| type == kFileName ) {
		fullPath = new char[ homeLen + len + 1 ];
		sprintf( fullPath, "%s%s%s", ( home ) ? home : "", ( home ) ? "/" : "", relPath );
	}
	else {
		fullPath = new char[ homeLen + len + 2 ];
		sprintf( fullPath, "%s%s%s/", ( home ) ? home : "", ( home ) ? "/" : "", relPath );
	}
	
	for ( int i = 0; i < strlen( fullPath ); ++i ) {
		if ( fullPath[i] == wrongSep ) {
			fullPath[i] = separator;
		}
	}
#	if defined (applec) || defined (powerc)
#		ifdef DEBUGPOWERC
	fprintf( stderr, "fullPath = %s\n", fullPath );
#		endif
#endif
	return fullPath;
}

Double Norm_mat( MatrixPtr ptr )
{
    int		i, k, rows = ptr->rows, cols = ptr->cols;
    Double	norm = 0.0;
    Double	**mat = ptr->mat;

    for ( k = 0; k < cols; ++k ) {
        for ( i = 0; i < rows; ++i ) {
            norm += MIN(fabs(mat[k][i]), 1.0e99) * MIN(fabs(mat[k][i]), 1.0e99);
        }
    }
    norm /= cols;
    norm = sqrt( norm );
    return norm;
}

Double MaxNorm_mat( MatrixPtr ptr, int *co, int *ro )
{
    int		i, k, rows = ptr->rows, cols = ptr->cols;
    Double	norm = 0.0;
	Double	**mat = ptr->mat;
	Double	absval;
	
	*co = 0;
	*ro = 0;

    for ( k = 0; k < cols; ++k ) {
        for ( i = 0; i < rows; ++i ) {
			absval = fabs( mat[k][i] );
			if ( absval > norm ) {
            	norm = absval;
				*co = k;
				*ro = i;
			}
        }
    }

    return norm;
}

#ifdef NEWMAXCHANGE
Double MaxChange( MatrixPtr oldPtr, MatrixPtr newPtr, int *co, int *ro )
{
    int		i, k, rows = oldPtr->rows, cols = oldPtr->cols;
    Double	norm = 0.0;
	Double	maxVec = 0.0;
	Double	**mat = oldPtr->mat;
	Double	**matNew = newPtr->mat;
	Double	absval = 0.0, fabmat, fabmatNew;

	if ( rows != newPtr->rows || cols != newPtr->cols ) {
		fprintf( stderr, "error in function 'MaxChange'\n" );
		exit( 2 );
	}

/* new procedure */
/*	int		maxLoc;*/
/*	Double	locVal;*/
/*	for ( i = 0; i < rows; ++i ) {*/
/*		maxVec = 0.0;*/
/*		for ( k = 0; k < cols; ++k ) {*/
/*			if ( ( locVal = fabs( mat[k][i] ) ) > maxVec )*/
/*				maxVec = locVal;*/
/*				maxLoc = k;*/
/*			}*/
/*		}*/
/*		fabmatNew = fabs( matNew[maxLoc][i] );*/
/*		if ( fabmatNew > 1.0e-8 && maxVec > 1.0e-8 && maxLoc / maxVec > 1e-5 ) {*/
/*			if ( fabmat > fabmatNew ) {*/
/*				absval = maxLoc / fabmatNew;*/
/*			}*/
/*			else {*/
/*				absval = fabmatNew / fabmat;*/
/*			}*/
/*		}*/
/*		if ( absval > norm ) {*/
/*			norm = absval;*/
/*			*co = k;*/
/*			*ro = i;*/
/*		}*/
/*    }*/
/* old procedure */
	for ( i = 0; i < rows; ++i ) {
		maxVec = 0.0;
		for ( k = 0; k < cols; ++k ) {
			maxVec = MAX( fabs( mat[k][i] ), maxVec );
		}
		if ( maxVec > 1.0e-08 ) {
			for ( k = 0; k < cols; ++k ) {
				fabmat = fabs( mat[k][i] );
				fabmatNew = fabs( matNew[k][i] );
				if ( fabmatNew > 1.0e-10 /*&& maxVec > 1.0e-5*/ && fabmat / maxVec > 9.0e-1 ) {
					if ( fabmat > fabmatNew ) {
						absval = fabmat / fabmatNew;
					}
					else {
						absval = fabmatNew / fabmat;
					}
				}
				if ( absval > norm ) {
					norm = absval;
					if ( norm > 1.0e5 ) {
						fprintf( stderr, "attention: %d\t%d\t%g\t%g\n", k, i, fabmat, fabmatNew );
						fprintf( stderr, "norm = %g\t%g\n", norm, fabmat / fabmatNew );
					}
					*co = k;
					*ro = i;
				}
			}
		}
	}

    return norm;
}
#else
Double MaxChange( MatrixPtr oldPtr, MatrixPtr newPtr, int *co, int *ro )
{
    int		i, k, rows = oldPtr->rows, cols = oldPtr->cols;
    Double	norm = 0.0;
	Double	maxVec = 0.0;
	Double	**mat = oldPtr->mat;
	Double	**matNew = newPtr->mat;
	Double	absval;
	MatrixPtr	locOlddeb = NewMatrix( oldPtr->phys_rows, oldPtr->phys_cols, kColumnPointers );
	copy_mat( locOlddeb, oldPtr );

	if ( rows != newPtr->rows || cols != newPtr->cols ) {
		fprintf( stderr, "error in function 'MaxChange'\n" );
		exit( 2 );
	}

	
	/*add_mat( oldPtr, newPtr );*/
    int	adrows = oldPtr->rows;
    int	adcols = oldPtr->cols;
    Double	val;
    for ( i = 0; i < adrows; ++i ) {
        for ( k = 0; k < adcols; ++k) {
			if ( fabs( matNew[k][i] ) >= fabs( mat[k][i] ) ) {
            	mat[k][i] += matNew[k][i];
				mat[k][i] *= 0.5;
			}
			else {
				val = 0.0;
				if ( matNew[k][i] ) {
					val += 1.0 / matNew[k][i];
				}
				if ( mat[k][i] ) {
					val += 1.0 / mat[k][i];
				}
            	mat[k][i] = ( val ) ? 1.0 / val : 0.0;
				mat[k][i] *= 2.0;
			}
		}
	}

//	mult_mat( oldPtr, 0.5 );
	div_mat( 2 , newPtr, oldPtr );
    for ( i = 0; i < rows; ++i ) {
		for ( k = 0; k < cols; ++k ) {
			absval = fabs( mat[k][i] );
			absval = ( absval < 1.0 && absval > 1.0e-2 ) ? 1.0 / absval : absval;
			if ( absval > norm && newPtr->mat[k][i] > 1.0e-60 && locOlddeb->mat[k][i] > 1.0e-60 ) {
//				fprintf( stderr, "fabs( mat[%d][%d] ) = %g\tabsval = %g\tY = %g\tYOld = %g\n", k, i
//							, fabs( mat[k][i] ), absval, newPtr->mat[k][i], locOlddeb->mat[k][i] );
				norm = absval;
				*co = k;
				*ro = i;
			}
        }
    }

//	if ( norm > 10 ) {
//		fprintf( stderr, "YN = %g\tval = %g\n", newPtr->mat[*co][*ro], fabs( mat[*co][*ro] ) );
//	}
/*	div_mat( 2 , newPtr, oldPtr );
    for ( i = 0; i < rows; ++i ) {
		maxVec = fabs( matNew[0][i] );
		for ( k = 1; k < cols; ++k ) {
			if ( fabs( matNew[k][i] ) > maxVec ) {
				maxVec = fabs( matNew[k][i] );
			}
		}
		for ( k = 0; k < cols; ++k ) {
			norm = ( fabs( matNew[k][i] ) > 0.01 * maxVec ) ?
						MAX( fabs( mat[k][i] ), norm ) : norm;
        }
    }*/

	DisposeMatrix( locOlddeb );
    return norm;
}
#endif
int MaxOfFabsVec( int len, Double *vec, int offset )
{
	int		maxPoint = 0;	
	int		off;
	Double	high;

	high = fabs( vec[0] );
	
	for ( int k = 1; k < len; ++k ) {
		if ( fabs( vec[off = k*offset] ) > high ) {
			high = fabs( vec[off] );
			maxPoint = k;
		}
	}
	
	return maxPoint;
}

Double ScaledNorm( MatrixPtr ptr, VectorPtr vecPtr )
{
    int		i, k, rows = ptr->rows, cols = ptr->cols;
    Double	norm = 0.0;
	Double	**mat = ptr->mat;
	Double	*scaler = vecPtr->vec;
	
	if ( rows > vecPtr->len ) {
		NewtonFatalError( "scaler not valid for current matrix" );
	}

    for ( k = 0; k < cols; ++k ) {
		for ( i = 0; i < rows; ++i ) {
            norm += mat[k][i] * mat[k][i] / ( scaler[i] * scaler [i] );
        }
    }
    norm /= cols;
    norm = sqrt( norm );
    return norm;
}

Double ScaledMaxNorm( MatrixPtr ptr, VectorPtr vecPtr )
{
	int		i, k, rows = ptr->rows, cols = ptr->cols;
    Double	norm = 0.0;
	Double	**mat = ptr->mat;
	Double	*scaler = vecPtr->vec;
	Double	locMax;
	
	if ( rows > vecPtr->len ) {
		NewtonFatalError( "scaler not valid for current matrix" );
	}

	for ( i = 0; i < rows; ++i ) {
		locMax = mat[0][i];
		for ( k = 1; k < cols; ++k ) {
            locMax = MAX( locMax, mat[k][i] );
        }
		norm += locMax / scaler[i];
    }
    return norm;
}

int copy_mat( MatrixPtr ptr1, MatrixPtr ptr2)
{
    int i, k, rows, cols;
	Double	**mat1 = ptr1->mat;
	Double	**mat2 = ptr2->mat;
    
    rows = ptr1->rows = ptr2->rows;
    cols = ptr1->cols = ptr2->cols;
	
	    
    for ( i = 0; i < rows; ++i ){
        for ( k = 0; k < cols; ++k){
            mat1[k][i] = mat2[k][i];
        }
    }

    return 0;
}

int copy_vec( VectorPtr ptr1, VectorPtr ptr2)
{
    int k;
    
	if ( ptr2->len > ptr1->phys_len ) {
		NewtonFatalError( "error: try to copy a vector to a shorter one" );
	}
	ptr1->len = ptr2->len;	
	    
    for ( k = 0; k < ptr2->len; ++k){
            ptr1->vec[k] = ptr2->vec[k];
    }

    return 0;
}

int add_mat( MatrixPtr ptr1, MatrixPtr ptr2 )
{
    int i, k, rows, cols;
	Double	**mat1 = ptr1->mat;
	Double	**mat2 = ptr2->mat;
    
    rows = ptr1->rows;
    cols = ptr1->cols;
    
    for ( i = 0; i < rows; ++i )
        for ( k = 0; k < cols; ++k)
            mat1[k][i] += mat2[k][i];
    return 0;
}

int sub_mat( int flag ,MatrixPtr ptr1, MatrixPtr ptr2 )
{
    /*  flag tells this function in which matrix solution should be stored      */
    /*  flag = 1: first matrix                                                  */
    /*  flag = 2: second matrix                                                 */
    
    int i, k, rows, cols;
    Double	**mat1 = ptr1->mat;
    Double	**mat2 = ptr2->mat;
    
    rows = ptr1->rows;
    cols = ptr1->cols;
    
    switch ( flag ){
        case 1:
            for ( i = 0; i < rows; ++i )
                for ( k = 0; k < cols; ++k)
                    mat1[k][i] -= mat2[k][i];
            break;
        case 2:
            for ( i = 0; i < rows; ++i )
                for ( k = 0; k < cols; ++k)
                    mat2[k][i] = mat1[k][i] - mat2[k][i];
            break;
        default: 
	    NewtonFatalError( "invalid call of function sub_mat" );    
    }
    return 0;
}

int div_mat( int flag ,MatrixPtr ptr1, MatrixPtr ptr2 )
{
    /*  flag tells this function in which matrix solution should be stored      */
    /*  flag = 1: first matrix                                                  */
    /*  flag = 2: second matrix                                                 */
    
    int i, k, rows, cols;
	Double	**mat1 = ptr1->mat;
	Double	**mat2 = ptr2->mat;
    
    rows = ptr1->rows;
    cols = ptr1->cols;
    
    switch ( flag ){
        case 1:
            for ( i = 0; i < rows; ++i ) {
                for ( k = 0; k < cols; ++k ) {
                    if ( mat2[k][i] ) {
						mat1[k][i]  /= mat2[k][i];
					}
					else {
						mat1[k][i] = 0.0;
					}
				}
			}
            break;
        case 2:
            for ( i = 0; i < rows; ++i ) {
                for ( k = 0; k < cols; ++k ) {
                    if ( mat2[k][i] ) {
						mat2[k][i] = mat1[k][i] / mat2[k][i];
					}
					else {
						mat2[k][i] = 0.0;
					}
				}
			}
            break;
        default: 
        	NewtonFatalError( "invalid call of function sub_mat" );    
    }
    return 0;
}

int mult_mat( MatrixPtr ptr, Double fact )
{
    int i, k, rows, cols;
	Double	**mat = ptr->mat;
	Double	*matk = NULL;
    
    rows = ptr->rows;
    cols = ptr->cols;
    
	for ( k = 0; k < cols; ++k) {
		matk = mat[k];
		for ( i = 0; i < rows; ++i ) {
			matk[i] *= fact;
		}
	}

    return 0;
}

void Test_convergence( TBVPSolverPtr solver, void *object )
{
	TNewtonPtr			bt = solver->bt;
	TTimePtr			time = solver->time;
	TDampPtr			damp = solver->damp;
	TContinuationPtr	        contin = solver->contin;
	MatrixPtr			locDy = bt->GetDy();
	if ( bt->GetTimedepFlag() ) {
		time->UnsetAdjustTime();
	}
	
	/*  compute norm.                                                           */
	bt->ScaleMatrix( locDy, bt->GetGrid()->GetSolutionScaler() );
	
	if ( bt->PrintFullRes() ) {
		FILE	*fpr = NULL;
		sprintf( bt->GetOutFileBuff(), "%sSolIncMat.tout", bt->GetOutputPath() );
		fpr = fopen( bt->GetOutFileBuff(), "w" );
		gPrnt->format = "%12.3E";
		gPrnt->title = "\nMatrix Scaled dy";
		fprintf( fpr, "niter = %d\n", bt->GetNIter() );
		PrintMatrix( locDy, gPrnt, fpr );
		fclose( fpr );
	}
		
	//printf("Norm of dY\n");
	bt->SetDyNorm( locDy );

	bt->UpdateRHS( locDy, object );
	
	bt->ScaleMatrix( locDy, bt->GetGrid()->GetSolutionScaler() );

        bt->SetResNorm( locDy );
	
	if ( bt->IsWriteResiduum() ) {
		bt->PrintResiduals( locDy, bt->GetVariableNames( object ) );
	}

	if ( bt->PrintFullRes() ) {
		FILE	*fpr = NULL;
		sprintf( bt->GetOutFileBuff(), "%sresidualmat.tout", bt->GetOutputPath() );
		if ( ( fpr = fopen( bt->GetOutFileBuff(), "w" ) ) == NULL ) {
			cerr << "#warning: can't open file residualmat.tout" << NEWL;
		}
		else {
/*			gPrnt->format = "%12.3E";
			gPrnt->title = "\nMatrix Scaled RHS";
			fprintf( fpr, "niter = %d\n", bt->GetNIter() );
			PrintMatrix( locDy, gPrnt, fpr );*/	
			Double	*x = bt->GetGrid()->GetCurrentGrid()->GetX()->vec;
			bt->PrintSolution( x, locDy->mat, bt->GetVariableNames( object ), fpr );
			fclose( fpr );
		}
	}

//  Test convergence of Newton iteration
	if ( bt->GetDyNorm() < bt->GetDyTol() || bt->GetResNorm() < bt->GetResTol() ){ // converged
		if ( bt->GetTimedepFlag() && bt->GetNIter() > 0 
			) {
			if ( time->GetTimeConverged() 
			&& (bt->GetGrid()->GetOneSolOneGrid()?(bt->GetGrid()->GetSolHasNewGrid()):TRUE)
				){ // final convergence
				bt->SetConvergeNewton();
			}
			else{
				bt->UnSetConvergeNewton();
				time->SetAdjustTime();
//				time->AdjustTimeStep( solver, TRUE, object ); // prepare for new time step
			}
		}
		else{
			bt->SetConvergeNewton();
		}
	}
	else{
		bt->UnSetConvergeNewton();
	}

	if ( bt->GetConvergNewton() == 1 || bt->GetNIter() >= bt->GetMaxIter() ||
		( ( !( bt->GetResNorm() < 1.0e10 ) || !( bt->GetDyNorm() < 1.0e10 ) ) 
			&& bt->GetNIter() > 3 ) || 
			( damp && ( damp->IsLambdaConvergent() == 0 ) ) ) {
//			( damp ) ? ( damp->IsLambdaConvergent() == 0 ) : NULL ) {
		bt->SetLeaveNewton();
	}
	else{
		bt->UnSetLeaveNewton();
	}
	
	if ( bt->GetTimedepFlag() && !time->GetTimeConverged()
			&& ( !( bt->GetResNorm() < 1.0e4 ) || !( bt->GetDyNorm() < 1.0e4 ) 
				|| bt->GetNIter() >= bt->GetMaxIter() ) ) {
		time->SetAdjustTime();
		bt->SetLeaveNewton();
//		time->AdjustTimeStep( solver, FALSE, object );
//		bt->UnSetLeaveNewton();
	}
	
	if ( bt->GetGrid()->GetOneSolOneGrid() && !bt->GetGrid()->GetSolHasNewGrid() && bt->GetConvergNewton() ) {
			bt->UnSetLeaveNewton();
	}

	if ( bt->GetContinuationFlag() && ( contin->GetParameter() < 1.0 
								|| bt->GetConvergNewton() == 0 ) ) {
		bt->UnSetLeaveContin();				
	}
	else{
		bt->SetLeaveContin();
	}

	if ( bt->GetNEquations() < bt->GetNVariables() && bt->GetConvergNewton() == TRUE ) {
		bt->UnSetLeaveContin();				
	}
}

void Monitor( TNewtonPtr bt, TContinuationPtr contin, TDampPtr damp, TTimePtr time )
{
	/*  Monitor.                                                        */
#ifdef NEWMONITOR
	if ( bt->GetContinuationFlag() ){
		fprintf( stderr, "ContPar = %-8g  ", contin->GetParameter() );
	}
	fprintf( stderr, "It = %-4d  GrdPts = %-4d  ",
	bt->GetNIter(), bt->GetGrid()->GetCurrentGrid()->GetNGridPoints() );
			
	if ( bt->GetDampFlag() ){
		fprintf( stderr, "Lambda = %-4.4g  ", damp->GetLambda() );
	} 
	if ( time && bt->GetTimedepFlag() && time->GetDeltaT() < time->GetDeltaTMax() ) {
		fprintf( stderr, "Delta_t = %-4g  ", time->GetDeltaT() );
	} 
	fprintf( stderr, "||dY|| = %-6g  ||Res|| = %-6g", bt->GetDyNorm(), bt->GetResNorm() );
	if ( bt->IsModified() ) {
		fprintf( stderr, "  modified\n" );
	}
	else {
		fprintf( stderr, "\n" );
	}
#else
	if ( bt->GetContinuationFlag() ){
		fprintf( stderr, "Cont. Param. = %-8g\t", contin->GetParameter() );
	}
	fprintf( stderr, "Iter = %-4d\tGridPts = %-4d\t",
	bt->GetNIter(), bt->GetGrid()->GetCurrentGrid()->GetNGridPoints() );
			
	if ( bt->GetDampFlag() ){
		fprintf( stderr, "Lambda = %-4.4g\t", damp->GetLambda() );
	} 
	if ( time && bt->GetTimedepFlag() && time->GetDeltaT() < time->GetDeltaTMax() ) {
		fprintf( stderr, "Delta_t = %-8g\t", time->GetDeltaT() );
	} 
	fprintf( stderr, "||dY|| = %-8g\t||Res|| = %-8g\t", bt->GetDyNorm(), bt->GetResNorm() );
	if ( bt->IsModified() ) {
		fprintf( stderr, "modified\n" );
	}
	else {
		fprintf( stderr, "\n" );
	}
#endif
}

void DoExit( TNewtonPtr bt, void *object )
{
	FILE	*fp;

	sprintf( bt->GetOutFileBuff(), "%sinterrupt%d", bt->GetOutputPath(), bt->GetNEquations() );
	fp = fopen( bt->GetOutFileBuff(), "w" );
	bt->WriteOutput( object, fp, "" );
	cerr << NEWL << "program stopped by user\noutput written to " << bt->GetOutFileBuff() << NEWL;
	fclose( fp );
	exit( 2 );
}

/*Double dfdyCentral( int i, int j, NodeInfoPtr nodeInfo )
{
    Double      fy;
	Double		*y = nodeInfo->y;
    Double      DeltaY = 0.0,
                DeltaY2 = 0.0;

    DeltaY = ( fabs( y[i] ) > 0.)
            ? .001 * y[i] : 1.e-8;
    DeltaY2 = 2. * DeltaY;

    y[i] += DeltaY;
    fy = flame->CalcFuncComponent( j, nodeInfo, TRUE );

    y[i] -= DeltaY2;
    fy -= flame->CalcFuncComponent( j, nodeInfo, TRUE );

    y[i] += DeltaY;
    fy /= DeltaY2;
    return fy;
}
*/

Double dfdyUpwind( int i, int j, CalcFuncComponentPtr func, NodeInfoPtr nodeInfo, void *object )
{	//	computes df_j / dy_i

    Double      fy;
	Double		*y = nodeInfo->y;
	Double 		storeY = y[i];
    Double      DeltaY = 0.0;
	
    DeltaY = ( fabs( y[i] ) > 0.)
            ? .00001 * y[i] : 1.e-14;

    DeltaY = y[i] += DeltaY;
    fy = func( j, nodeInfo, object, TRUE );

    y[i] = storeY;
    fy -= func( j, nodeInfo, object, kDoNothing );

    fy /= ( DeltaY - storeY );
	
    return fy;
}

/*Double dfdyUpwind( int i, int j, NodeInfoPtr nodeInfo )
{	//	computes df_j / dy_i

    Double      fy;
	Double		*y = nodeInfo->y;
	Double 		storeY = y[i];
    Double      DeltaY = 0.0;
	
    DeltaY = ( fabs( y[i] ) > 0.)
            ? .00001 * y[i] : 1.e-14;

    DeltaY = y[i] += DeltaY;
    fy = flame->CalcFuncComponent( j, nodeInfo, TRUE );

    y[i] = storeY;
    fy -= flame->CalcFuncComponent( j, nodeInfo );

    fy /= ( DeltaY - storeY );
	
    return fy;
}
*/

Double Integrate( Double *x, Double *y, Double x1, Double x2, int nPoints, int &i1, int &i2 )
{
	int i;
	for ( i = 0; i < nPoints; ++i ) {
		if ( x1 <= x[i] ) {
			i1 = i;
			break;
		}
	}
	for ( ; i < nPoints; ++i ) {
		if ( x2 <= x[i+1] ) {
			i2 = i;
			break;
		}
	}
	
	return Integrate( x, y, i1, i2 );
}

Double Integrate( Double *x, Double *y, int i1, int i2 )
{
	Double	integral = y[i1] * ( x[i1+1] - x[i1] ) + y[i2] * ( x[i2] - x[i2-1] );
	
	for ( int i = i1+1; i < i2-1; ++i ) {
		integral += y[i] * ( x[i+1] - x[i-1] );
	}
	
	return 0.5 * integral;
}

void NewtonFatalError(const char *str)
{
	fprintf(stderr, "### Error: %s.\n", str);
	exit(2);
}

void CopyMatrix( MatrixPtr source, MatrixPtr dest )
{
	//  make a copy of the physical domain and adjust the logical domain

	if (( source->phys_rows == dest->phys_rows ) && 
		( source->phys_cols == dest->phys_cols ) &&
		( source->partition == dest->partition )) {
		memcpy( dest->mat[0], source->mat[0], source->phys_rows * source->phys_cols * sizeof( Double ));
		dest->rows = source->rows;
		dest->cols = source->cols;
	}
	else {
		cerr << "#error: CopyMatrix failed, source and destination have different "
			 << "physical size or partition type" << NEWL;
		exit( 2 );
	}
}
