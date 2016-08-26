#include "Newton.h"
#include "Spline.h"

#undef SPLINE
#undef SAVEM
#undef SAVEF
#undef DEBUGGRID
#undef DEBUGSCALER
#undef DEBUGESTIMATEDERROR
#undef DEBUGSCALERANDY
#undef NEWWEIGHTS

//GB
#undef SLOWGRID

enum Smooth { kMax3, kMax5, kSmooth3, kSmooth5 };

void TAdaptiveGrid::new_grid( TBVPSolverPtr solver, void *object )
{
    register int    i, m;
             int    n = fine->GetNGridPoints(), 
					n_old = n;
    register Double f, a, p = q;
             Double *localX = fine->GetX()->vec, 
                    *x2 = aux->vec, 
                    *locM = M->vec;
	TNewtonPtr		bt = solver->bt;
	TDampPtr		damp = solver->damp;
	TTimePtr		time = solver->time;
	int				maxGridPoints = bt->GetMaxGridPoints();
    
#ifdef HP
	if ( !bt->IsGriddingStep() && gExit != kNewGrid ) {
		return;
	}
#else
	if ( !bt->IsGriddingStep() ) {
		return;
	}
#endif
	else {
#ifdef HP
		if ( gExit == kNewGrid ) {
			InstallInterruptHP();
		}
#endif
	}
	fprintf( stderr, "generate new grid\n");

/*  solution on coarse grid  */
	Make_coarse_grid( bt );
	switch_to_coarse( solver );
	bt->fUtFuncs->PostIter( object );

#ifdef DEBUGGRID
	FILE *fpCoarse = fopen( "coarsegrid.dout", "w" );
	bt->PrintSolution( coarse->GetX()->vec, coarse->GetY()->mat, bt->GetVariableNames( object ), fpCoarse );
	fclose( fpCoarse );
#endif
#ifdef NEWWEIGHTS
	CopyMatrix( bt->GetDy(), bt->GetDySaved() );	// save rhs of the coarse grid
#else
	bt->Newton_step( damp, time, object );
#endif

#ifdef DEBUGGRID
	FILE *fpCoarseSolu = fopen( "coarsesolution.dout", "w" );
	bt->PrintSolution( coarse->GetX()->vec, coarse->GetY()->mat, bt->GetVariableNames( object ), fpCoarseSolu );
	fclose( fpCoarseSolu );
#endif
	switch_to_fine( solver );
#ifdef NEWWEIGHTS
	bt->fUtFuncs->PostIter( object );
#endif

	compute_M( bt );
	compute_F( n_old, localX, locM, bt );
	if ( bt->GetInitialGridpoints() < 0 ) {
		n = -bt->GetInitialGridpoints();
		p = locM[n_old] / ( ( Double )n + 1.0);       /* p = F(a_N) / (N+1) */
	} 
	else if (p < 0.0) {                                  /* don't change the number of nodes */
		n = bt->GetMaxGridPoints();
		p = locM[n_old] / ( ( Double )n + 1.0);       /* p = F(a_N) / (N+1) */
	} 
	else {
	        n = (int) MIN( maxGridPoints, ( locM[n_old] / p - 1 ) );
		n = ( n + n_old ) / 2;

		if ( !( n % 2 ) ) ++n;
                                                    /* odd grid points for fine grid required */
		if ( n > maxGridPoints )
		  n = maxGridPoints;                 /* too much accuracy requested, use n_max points */
		p = locM[n_old] / ((Double)n + 1.0);           /* correct p */
	}
	for ( i = 0, m = 0; i < n; ++i ) {
	  f = (Double)(i+1) * p;                          
	  if ( f <= locM[0] ){
            a = f / locM[0];
            x2[i] = bt->GetLeft() +  a * (localX[0] - bt->GetLeft());
#ifdef DEBUGGRID
			//		fprintf(stdout, "i = %d\tm = %d\t%g\t%g\t%g\n", i, 0,  locM[m], f, locM[m+1]);
#endif
		}
        else if ( f >= locM[n_old-1] ){
            a = (f - locM[n_old-1]) / (locM[n_old] - locM[n_old-1]);
            x2[i] = localX[n_old-1] +  a * (bt->GetRight() - localX[n_old-1]);
#ifdef DEBUGGRID
			//		fprintf(stdout, "i = %d\tm = %d\t%g\t%g\t%g\n", i, n_old-1,  locM[m], f, locM[m+1]);
#endif
		}
		else{
            while (f >= locM[++m]){
                if ( m > n_old ){
                    fprintf( stderr, "#error during grid generation\nf = %g\ti=%d\tm = %d\n", f, i, m);
                    exit(2);
                }
            }
            --m;
#ifdef DEBUGGRID
			//		fprintf(stdout, "i = %d\tm = %d\t%g\t%g\t%g\n", i, m,  locM[m], f, locM[m+1]);
#endif
            a = (f - locM[m]) / (locM[m+1] - locM[m]);
            x2[i] = localX[m] +  a * (localX[m+1] - localX[m]);
        }
    }
/*  store y_old in auxmat
    store x_old in locM */
    for ( i = 0; i < n_old; ++i ){
	   locM[i] = localX[i];
        for ( m = 0; m < bt->GetNVariables(); ++m ){
            auxmat->mat[i][m] = fine->GetY()->mat[i][m];
        }
    }

#ifdef DEBUGGRID
	FILE *fp = fopen( "vorInterpolation", "w" );
	bt->PrintSolution( localX, fine->GetY()->mat, bt->GetVariableNames( object ), fp );
	fclose( fp );
#endif

// set new grid ( n, x, y )
#ifndef SLOWGRID
	for (i = 0; i < n; ++i) localX[i] = x2[i];
#else
	Double ratio = 0.5;
	for (i = 0; i < n; ++i)
	  localX[i] = (1.0-ratio) * locM[i] + ratio * x2[i];
#endif
	
	fine->AdjustNGridPoints( n );
	
	smooth_x( localX, aux->vec, fine->GetNGridPoints(), bt );       /* smooth the new grid */
	
#ifdef SPLINE
	SplineInterpolateY( bt, locM, auxmat->mat, n_old, localX, fine->GetY()->mat );
#else
	LinearInterpolateY( bt, locM, auxmat->mat, n_old, localX, fine->GetY()->mat );
#endif
	
	if ( bt->GetDampFlag() ) {
		copy_mat( auxmat, damp->GetDyOld() );
#ifdef SPLINE
		SplineInterpolateY( bt, locM, auxmat->mat, n_old, localX, damp->GetDyOld()->mat );
#else
		LinearInterpolateY( bt, locM, auxmat->mat, n_old, localX, damp->GetDyOld()->mat );
#endif
	}
	if ( bt->GetTimedepFlag() ) {
		copy_mat( auxmat, time->GetYOld() );
		LinearInterpolateMat( bt, locM, auxmat->mat, n_old, localX, time->GetYOld()->mat );

		copy_mat( auxmat, time->GetYPrime() );
		LinearInterpolateMat( bt, locM, auxmat->mat, n_old, localX, time->GetYPrime()->mat );

		copy_vec( time->GetXOld(), fine->GetX() );
	}

	solver->UpdateAllDimensions( n );
	
	bt->fUtFuncs->PostIter( object );

	fSolHasNewGrid = TRUE;
	fIterOfGeneration = bt->GetNIter();
	
#ifdef DEBUGGRID
	TContinuationPtr contin = NULL;
	
	Test_convergence( solver, object );
	Monitor( bt, contin, damp, time );

	fp = fopen( "nachInterpolation", "w" );
	bt->PrintSolution( localX, fine->GetY()->mat, bt->GetVariableNames( object ), fp );
	fclose( fp );
#endif
}

void TAdaptiveGrid::ToOldSolution( TBVPSolverPtr solver, VectorPtr xOld, MatrixPtr yOld )
{
	int gridPoints = xOld->len;

	copy_vec( fine->GetX(), xOld );
	copy_mat( fine->GetY(), yOld );

    fine->AdjustNGridPoints( gridPoints );
	solver->UpdateAllDimensions( gridPoints );
}

void TAdaptiveGrid::Make_coarse_grid( TNewtonPtr bt )
{
	//	|--*--*--*--|	fine grid
	//	|-----*-----|	caorse grid

	int	coarseGridPoints = ( fine->GetNGridPoints() - 1 ) / 2;

    if ( !(bt->GetMaxGridPoints() % 2) ){
        fprintf( stderr, "#error while generating coarse grid: number of nodes of fine grid has to be odd\n" );
        exit(2);
    }

	coarse->AdjustNGridPoints( coarseGridPoints );

	//	Set up local variables for easier access to arrays and matrices
	Double	*x_coarse = coarse->GetX()->vec;	// pointer to x-values of coarse grid
	Double	*x_fine   = fine->GetX()->vec;		// pointer to x-values of fine grid
	Double	**m_coarse = coarse->GetY()->mat;	// dependent variables on coarse grid
	Double	**m_fine   = fine->GetY()->mat;		// dependent variables on fine grid
	Double	*yLeftFine = fine->GetYLeft()->vec;
	Double	*yLeftCoarse = coarse->GetYLeft()->vec;
	Double	*yRightFine = fine->GetYRight()->vec;
	Double	*yRightCoarse = coarse->GetYRight()->vec;
	int		i;							// inner loop counter (function scope!)
	int		variables = bt->GetNVariables();	
	
    for ( int k = 0; k < coarse->GetNGridPoints(); ++k ) {
		x_coarse[k] = x_fine[2*k+1];
        for ( i = 0; i < variables; ++i ) {
		   m_coarse[k][i] = m_fine[2*k+1][i];
        }
    }
	//	copy boundary conditions of fine mesh
	copy( variables, fine->GetBcLeft()->vec, 1, coarse->GetBcLeft()->vec, 1 );
	copy( variables, fine->GetBcRight()->vec, 1, coarse->GetBcRight()->vec, 1 );
	copy( variables, fine->GetYLeft()->vec, 1, coarse->GetYLeft()->vec, 1 );
	copy( variables, fine->GetYRight()->vec, 1, coarse->GetYRight()->vec, 1 );
/*    for ( i = 0; i < variables; ++i ) {
		yLeftCoarse[i] = yLeftFine[i];
		yRightCoarse[i] = yRightFine[i];
	}*/
}

void TAdaptiveGrid::UpdateGridDimension( void )
{
	int	number = fine->GetNGridPoints();
	
	auxmat->cols = number;
	aux->len = number;
	aux2->len = number;
	M->len = number+1;

	fine->UpdateOneGridDimension( number );
	number = ( fine->GetNGridPoints() - 1 ) / 2;
	coarse->UpdateOneGridDimension( number );
}

void TAdaptiveGrid::SplineInterpolateY( TNewtonPtr bt, Double *x_old, Double **y_old, int n_old,
						Double *x, Double **y )
{
	int			i, k;
	TGridPtr	theGrid =  bt->GetGrid()->GetCurrentGrid();
	SplinePtr	theSpline = NULL;
	int			nGridPoints = theGrid->GetNGridPoints();
	int			nEquations = bt->GetNVariables();
	Double  	*yLeft = theGrid->GetYLeft()->vec;
	Double		*yRight = theGrid->GetYRight()->vec;
	Double		*xVec = aux->vec;	// xVec[0] contains left boundary
	Double		*yVec = aux2->vec;	// yVec[0] contains left BC
	Double		left = bt->GetLeft();
	Double		right = bt->GetRight();
	Double		leftSlope;
	Double		rightSlope;
	
	xVec[0] = left;
	for ( k = 0; k < n_old; ++k ) {
		xVec[k+1] = x_old[k];
	}
	xVec[n_old+1] = right;
	
	for ( i = 0; i < nEquations; ++i ) {	// loop variables
		yVec[0] = yLeft[i];
		for ( k = 0; k < n_old; ++k ) {		// copy current variable to workspace
			yVec[k+1] = y_old[k][i];
			
		}
		yVec[n_old+1] = yRight[i];
		
		leftSlope = ( yVec[1] - yVec[0] ) / ( xVec[1] - xVec[0] );
		rightSlope = ( yVec[n_old+1] - yVec[n_old] ) / (  xVec[n_old+1] - xVec[n_old] );
		theSpline = ComputeSimpleSpline( xVec, yVec, n_old+2, FALSE, leftSlope, FALSE, rightSlope, NULL, TRUE );
		SplineInterpolate( theSpline, x, yVec, nGridPoints );
		FreeSpline( theSpline );

		for ( k = 0; k < nGridPoints; ++k ) {	// copy workspace to vector of solution
			y[k][i] = yVec[k];
		}
	}
//	FreeSpline( theSpline );
}

void TAdaptiveGrid::LinearInterpolateY( TNewtonPtr bt, Double *x_old, Double **y_old, int n_old,
						Double *x, Double **y )
{
    int     	k, i, m;
	TGridPtr	theGrid =  bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = theGrid->GetNGridPoints();
	int			nEquations = bt->GetNVariables();
	Double  	*yleft = theGrid->GetYLeft()->vec;
	Double		*yright = theGrid->GetYRight()->vec;

    for ( k = 0, m = 0; k < nGridPoints; ++k ) {
        if ( x[k] <= x_old[0] ) {
            for (i = 0; i < nEquations; ++i ) {
                y[k][i] = yleft[i] +
                    ( y_old[0][i] - yleft[i] ) *
                    ( x[k] - bt->GetLeft() ) /
                    ( x_old[0] - bt->GetLeft() );
            }
        }
        else if ( x[k] >= x_old[n_old-1] ){
            for (i = 0; i < nEquations; ++i ){
                y[k][i] = y_old[n_old-1][i] +
                    ( yright[i] - y_old[n_old-1][i] ) *
                    ( x[k] - x_old[n_old-1] ) /
                    ( bt->GetRight() - x_old[n_old-1] );
            }
        }
        else{
            while( x_old[++m] <= x[k] ){
				if ( m > n_old+1 ){
					fprintf( stderr, "#error while interpolating y on new grid\n");
					exit(2);
				}
			}
			--m;
            for ( i = 0; i < nEquations; ++i ){
                y[k][i] = y_old[m][i] +
                    ( y_old[m+1][i] - y_old[m][i] ) *
                    ( x[k] - x_old[m] ) /
                    ( x_old[m+1] - x_old[m] );
            }
        }
    }
}

void TAdaptiveGrid::LinearInterpolateMat( TNewtonPtr bt, Double *x_old, Double **y_old, int n_old,
						Double *x, Double **y )
{
    int     	k, i, m;
	TGridPtr	theGrid =  bt->GetGrid()->GetCurrentGrid();
	int			nGridPoints = theGrid->GetNGridPoints();
	int			nEquations = bt->GetNVariables();

    for ( k = 0, m = 0; k < nGridPoints; ++k ) {
        if ( x[k] <= x_old[0] ) {
            for (i = 0; i < nEquations; ++i ) {
                y[k][i] = y_old[0][i] +
                    ( y_old[1][i] - y_old[0][i] ) *
                    ( x[k] - x[0] ) /
                    ( x_old[1] - x_old[0] );
            }
        }
        else if ( x[k] >= x_old[n_old-1] ){
            for (i = 0; i < nEquations; ++i ){
                y[k][i] = y_old[n_old-1][i] +
                    ( y_old[n_old-1][i] - y_old[n_old-2][i] ) *
                    ( x[k] - x_old[n_old-1] ) /
                    ( x_old[n_old-1] - x_old[n_old-2] );
            }
        }
        else{
            while( x_old[++m] <= x[k] ){
				if ( m > n_old+1 ){
					fprintf( stderr, "#error while interpolating y on new grid\n");
					exit(2);
				}
			}
			--m;
            for ( i = 0; i < nEquations; ++i ){
                y[k][i] = y_old[m][i] +
                    ( y_old[m+1][i] - y_old[m][i] ) *
                    ( x[k] - x_old[m] ) /
                    ( x_old[m+1] - x_old[m] );
            }
        }
    }
}

void TAdaptiveGrid::compute_M( TNewtonPtr bt )
{
//	estimates truncation error

	register int    i, j, 
	   				n = fine->GetNGridPoints(), 
					n_coarse = coarse->GetNGridPoints();
	int				nEquation = bt->GetNEquations();
	Double          *locM = M->vec, 
	                *x = fine->GetX()->vec,
#ifdef NEWWEIGHTS
	                **yFine = bt->GetDy()->mat,			// use residuals instead of solution
	                **yCoarse = bt->GetDySaved()->mat,
#else
	                **yFine = fine->GetY()->mat,
	                **yCoarse = coarse->GetY()->mat,
#endif
	                **error = auxmat->mat;

#ifdef SAVEM
    static int  m = 0;
    char        fname[32];
    FILE        *fp;
#endif

//	ScaleSolution( fine );
//	ScaleSolution( coarse );

#ifdef DEBUGSCALERANDY
    FILE        *fpy;
    char        fgame[32];

    sprintf( fgame, "YScaled%d", bt->GetNIter() );
    fpy = fopen( fgame, "w" );
    gPrnt->title = "\nVector scaler";
    PrintVector( fSolutionScaler, gPrnt, fpy );
    gPrnt->title = "\nMatrix y_scaled";
    PrintMatrix( grid->GetY(), gPrnt, fpy );
    fclose(fpy);
#endif

#ifdef DEBUGESTIMATEDERROR
 	{   
		FILE        *fpe;
    	char        fename[32];

	    sprintf( fename, "truncerror%d.dout", bt->GetNIter() );
	    fpe = fopen( fename, "w" );
		
	    int         i, k;
	
	    fprintf( fpe, "*\n" );
		
		fprintf( fpe, "%-12s", "eta" );
		for ( i = 0; i < nEquation; ++i ) {
			fprintf( fpe, "\t%-12s", "var" );
		}

	    for ( k = 0; k < n_coarse; ++k ){
	        fprintf( fpe, "\n%-12E", x[2*k+1] );
	        for ( i = 0; i < nEquation; ++i ){
	            fprintf( fpe, "\t%-12E", fabs( yFine[2*k+1][i] - yCoarse[k][i] ) );
	        }
	    }
	    fprintf( fpe, "\n" );

	    fclose(fpe);
	}
#endif

// compute error
    for (i = 0; i < n_coarse; ++i) {
		for ( j = 0; j < nEquation; ++j ){
			error[2*i+1][j] = fabs( yFine[2*i+1][j] - yCoarse[i][j] );
/*			if ( j != 86 && j != 87 ) {*/
/*				error[2*i+1][j] = fabs( yFine[2*i+1][j] - yCoarse[i][j] );*/
/*			}*/
/*			else {*/
/*				error[2*i+1][j] = 0.0;*/
/*			}*/
		}
		locM[2*i+1] /= nEquation;
	}

// scale it
	for ( i = 0; i < nEquation; ++i ) {
		integral[i] = error[1][i] * ( x[1] - bt->GetLeft() );
		for ( j = 1; j < n_coarse; ++j ) {
			integral[i] += ( error[2*j+1][i] + error[2*j-1][i] ) * ( x[2*j+1] - x[2*j-1] );
		}
		integral[i] += error[2*n_coarse - 1][i] * (  bt->GetRight() - x[n_coarse-1] );
		integral[i] *= 0.5;
		if ( integral[i] / ( bt->GetRight() - bt->GetLeft() ) < 1.0e-20 ) {
			integral[i] = 1.0;
		}
	}

#ifdef DEBUGESTIMATEDERROR
 	{   
		FILE        *fpe;
    	char        fename[32];

	    sprintf( fename, "truncerrorscaled%d.dout", bt->GetNIter() );
	    fpe = fopen( fename, "w" );
		
	    int         i, k;
	
	    fprintf( fpe, "*\n" );
		
		fprintf( fpe, "%-12s", "eta" );
		for ( i = 0; i < nEquation; ++i ) {
			fprintf( fpe, "\t%-12s", "var" );
		}

	    for ( k = 0; k < n_coarse; ++k ){
	        fprintf( fpe, "\n%-12E", x[2*k+1] );
	        for ( i = 0; i < nEquation; ++i ){
	            fprintf( fpe, "\t%-12E", fabs( yFine[2*k+1][i] - yCoarse[k][i] ) / integral[i] );
	        }
	    }
	    fprintf( fpe, "\n" );

	    fclose(fpe);
	}
#endif

// error = sum_i{ ( y_iFine - y_iCoarse ) / int( y_iFine - y_iCoarse ) }
	for (i = 0; i < n_coarse; ++i) {
		locM[2*i+1] = 0.0;
		for ( j = 0; j < nEquation; ++j ){
			locM[2*i+1] += error[2*i+1][j] / integral[j];
		}
		locM[2*i+1] /= nEquation;
	}
	locM[0] = ( x[0] - bt->GetLeft() ) * locM[1];
	locM[0] /= x[1] - bt->GetLeft();
	locM[n-1] = ( bt->GetRight() - x[n-1] ) * locM[n-2];
	locM[n-1] /= bt->GetRight() - x[n-2];

// linear interpolation of error to even gridpoints
	for (i = 2; i < n-2; i += 2) {
	   locM[i] = ( x[i+1] - x[i] ) * locM[i-1] + ( x[i] - x[i-1] ) * locM[i+1];
	   locM[i] /= x[i+1] - x[i-1];
	}
    
	/* filter the raw error estimates */
	filter( kMax5 );
	//    filter( kSmooth3 );

#ifdef SAVEM
    sprintf( fname, "M%02d.dout", m++ );
    fp = fopen( fname, "w" );
    fprintf( fp, "*\n" );
    fprintf( fp, "x\tM\n" );
    for (i = 0; i < n; ++i) fprintf( fp, "%g\t%g\n", x[i], locM[i] );
    fclose(fp);
#endif
	
//	find minimum
	if ( fAlpha ) {
		int		i1, i2; // i1 and i2 are set in function Integral
		Double	dM;
	
		dM = ( fAlpha * Integrate( x, locM, 0, n ) 
				- Integrate( x, locM, fStart, fEnd, n, i1, i2 ) ) 
				/ ( ( 1.0 - fAlpha ) * ( fEnd - fStart ) );
				
/*		cerr << "Integrate( x, locM, 0, n )M = " << Integrate( x, locM, 0, n ) << NEWL;
		cerr << "Integrate( x, locM, fStart, fEnd, n, i1, i2 ) = " << Integrate( x, locM, fStart, fEnd, n, i1, i2 ) << NEWL;
		cerr << "dM = " << dM << NEWL;
		cerr << "i1 = " << i1 << NEWL;
		cerr << "i2 = " << i2 << NEWL;*/
			
		for ( i = i1; i < i2; ++i ) {
			locM[i] += dM;
		}
	}
    filter( kSmooth3 );

	
#ifdef SAVEM
    sprintf( fname, "Mnew%02d.dout", m-1 );
    fp = fopen( fname, "w" );
	fprintf( fp, "*\n" );
    fprintf( fp, "x\tM\n" );
    for (i = 0; i < n; ++i) fprintf( fp, "%g\t%g\n", x[i], locM[i] );
    fclose(fp);
#endif

//	DeScaleSolution( fine );
//	DeScaleSolution( coarse );
}

void TAdaptiveGrid::DeScaleSolution( TGridPtr theGrid, Double **y )
{
	int		i, k;
	int		gridPoints = theGrid->GetNGridPoints();
	int		variables = theGrid->GetY()->rows;
//	Double	**y = theGrid->GetY()->mat;
	Double	*scaler = fSolutionScaler->vec;

	for ( i = 0; i < variables; ++i ) {
		for ( k = 0; k < gridPoints; ++k ) {
			y[k][i] *= scaler[i];
		}
	}
}

void TAdaptiveGrid::ScaleSolution( TGridPtr theGrid, Double **y )
{
	int		i, k;
	int		gridPoints = theGrid->GetNGridPoints();
	int		variables = theGrid->GetY()->rows;
//	Double	**y = theGrid->GetY()->mat;
	Double	*scaler = fSolutionScaler->vec;

	for ( i = 0; i < variables; ++i ) {
		for ( k = 0; k < gridPoints; ++k ) {
			y[k][i] /= scaler[i];
		}
	}
}

/*void TAdaptiveGrid::SetSolutionScaler( void )
{
	int		i, k;
	int		gridPoints = fine->GetNGridPoints();
	int		variables = fine->GetY()->rows;
	Double	**y = fine->GetY()->mat;
	Double	*x = fine->GetX()->vec;
	Double	*yLeft = fine->GetYLeft()->vec;
	Double	*yRight = fine->GetYRight()->vec;
	Double	*scaler = fSolutionScaler->vec;
	
	for ( i = 0; i < variables; ++i ) {
		scaler[i] = ( yLeft[i] + y[0][i] ) * ( x[0] - GetLeft() );
		for ( k = 1; k < gridPoints; ++k ) {
			scaler[i] += ( y[k][i] + y[k-1][i] ) * ( x[k] -  x[k-1] );
		}
		scaler[i] += ( y[gridPoints-1][i] + yRight[i] ) * ( GetRight() - x[gridPoints-1] );
		scaler[i] *= 0.5;
	}
#ifdef DEBUGSCALER
	FILE        *fpy;

	fpy = fopen( "scaler.tout", "w" );
	PrintVector( fSolutionScaler, gPrnt, fpy );
	fclose(fpy);
#endif
}
*/

void TAdaptiveGrid::SetSolutionScaler( void )
{
	int		i, k;
	int		gridPoints = fine->GetNGridPoints();
	int		variables = fine->GetY()->rows;
	Double	**y = fine->GetY()->mat;
	Double	*yLeft = fine->GetYLeft()->vec;
	Double	*yRight = fine->GetYRight()->vec;
	Double	*scaler = fSolutionScaler->vec;
	Double	locMax;
	
	for ( i = 0; i < variables; ++i ) {
		for ( k = 0, locMax = fabs( yLeft[i] ); k < gridPoints; ++k ) {
			locMax = MAX( locMax, fabs( y[k][i] ) );
		}
		scaler[i] = MAX( fabs( yRight[i] ), locMax );
		if ( scaler[i] < 1.0e-10 ) {
			scaler[i] = 1.0;
		}
	}

#ifdef DEBUGSCALER
	FILE        *fpy;

	fpy = fopen( "scaler.tout", "w" );
	PrintVector( fSolutionScaler, gPrnt, fpy );
	fclose(fpy);
#endif
}

void TAdaptiveGrid::compute_F( int n, Double *x, Double *theM, TNewtonPtr bt )
{
    /*  compute F(a_i) and store it in M(i), where F_left is 0.0 
        and F_right is stored in M[n]  */

    register int    i;
    register Double c, w, w_old, Mmax, Mmin;

#ifdef SAVEF
    static int  m = 0;
    char        fname[32];
    FILE        *fp;
#endif
    
    Mmin = Mmax = theM[0];             /* get min and max of M */
    for ( i = 1; i < n; ++i ) {
        if (theM[i] > Mmax)    Mmax = theM[i];
        if (theM[i] < Mmin)    Mmin = theM[i];
    }
    if (R >= Mmax / Mmin) R = 0.999 * Mmax / Mmin;
    c = (R - 1.0) / (Mmax - Mmin * R);
    w = w_old = 1.0 + c * theM[0];             /* 1st weight */
    theM[0] = 0.5 * ( w_old + w ) * ( x[0] - bt->GetLeft() );
    for ( i = 1; i < n; ++i ) {
        w_old = w;
        w = 1.0 + c * theM[i];
        theM[i] = theM[i-1] + 0.5 * (w_old + w) * (x[i] - x[i-1]);
    }
    w = 1.0 + c * theM[n];
    theM[n] = theM[n-1] + w_old * (bt->GetRight() - x[n-1]);
    
    
#ifdef SAVEF
    sprintf( fname, "F%02d.dat", m++ );
    fp = fopen( fname, "w" );
    fprintf( fp, "x\tF\n" );
    for (i = 0; i < n; ++i) fprintf( fp, "%g\t%g\n", x[i], theM[i] );
	fprintf( fp, "1.0\t%g\n", theM[n] );
    fclose(fp);
#endif
}

void TAdaptiveGrid::max3Filter( Double *theM, Double *auxarray, int n ) /* use a symmetric 3 point stencil */
{
    register int    i, nm1 = n - 1;
    register Double current_max;
    
    memcpy( auxarray, theM, n * sizeof(Double) );

    theM[0] = MAX( auxarray[0], auxarray[1] );
    for ( i = 1; i < nm1; ++i ) {
        current_max = MAX(auxarray[i-1], auxarray[i]);
        theM[i] = MAX(current_max, auxarray[i+1]);
    }
    theM[nm1] = MAX(auxarray[nm1-1], auxarray[nm1]);
}


void TAdaptiveGrid::max5Filter( Double *x, Double *auxarray, int n )
{
    register int i, n2 = n - 2;
    register Double max1, max2;
    
    memcpy( auxarray, x, n * sizeof(Double) );
    
    max1 = MAX( auxarray[1], auxarray[2] );
    x[0] = MAX( auxarray[0], max1 );
    max1 = MAX( auxarray[0], auxarray[1] );
    max2 = MAX( auxarray[2], auxarray[3] );
    x[1] = MAX( max1, max2 );
    for ( i = 2; i < n2; ++i ) {
        max1 = MAX( auxarray[i-2], auxarray[i-1] );
        max2 = MAX( auxarray[i+1], auxarray[i+2] );
        max1 = MAX( max1, max2 );
        x[i] = MAX( auxarray[i], max1 );
    }
    max1 = MAX( auxarray[n-3], auxarray[n-2] );
    x[n-1] = MAX( auxarray[n-1], max1 );
    max1 = MAX( auxarray[n-4], auxarray[n-3] );
    max2 = MAX( auxarray[n-1], auxarray[n-2] );
    x[n-2] = MAX( max1, max2 );

}


void TAdaptiveGrid::smoothFilter3( Double *x, Double *auxarray, int n )
{
    register int i, m = n-1;
    register Double aThird = 1.0 / 3.0;
    
    memcpy( auxarray, x, n * sizeof(Double) );
    
    x[0] = 0.5 * (auxarray[0] + auxarray[1]);
    for ( i = 1; i < m; ++i )
        x[i] = (auxarray[i-1] + auxarray[i] + auxarray[i+1]) * aThird;
    x[m] = 0.5 * (auxarray[m-1] + auxarray[m]);
}


void TAdaptiveGrid::filter( int which )
{
    /*  The switch statement is clumsy and will eventually be
        replaced by a lookup table. */
    switch ( which ) {
        case kMax3:
            max3Filter( M->vec, aux->vec, fine->GetNGridPoints() );
            break;
        case kMax5:
            max5Filter( M->vec, aux->vec, fine->GetNGridPoints() );
            break;
        case kSmooth3:
            smoothFilter3( M->vec, aux->vec, fine->GetNGridPoints() );
            break;
        default:
            fprintf( stderr, "### wrong filter option: ID = %d\n", which );
            exit ( 2 );
    }
}

void TAdaptiveGrid::smooth_x( Double *x, Double *auxarray, int n, TNewtonPtr bt )    /* weakly smooth independant variable */
{
    register int i;
    register Double weight = 0.25;
    
    memcpy( auxarray, x, n * sizeof(Double) );

    x[0] = ( bt->GetLeft() + 2.0 * auxarray[0] + auxarray[1] ) * weight;
    x[n-1] = ( auxarray[n-2] + 2.0 * auxarray[n-1] + bt->GetRight() ) * weight;
    for ( i = 1; i < n-1; ++i )
        x[i] = ( auxarray[i-1] + 2.0 * auxarray[i] + auxarray[i+1] ) * weight;
}

void TAdaptiveGrid::InitGridding( TNewtonPtr bt, TBVPInputPtr input )
{
   	int	maxGridPoints = bt->GetMaxGridPoints();
	int variables = bt->GetNVariables();

	R = input->fR;
	q = input->fQ;
	fStart = input->fStart;
	fEnd = input->fEnd;
	fAlpha = input->fAlpha;
	
	fOneSolOneGrid = input->fOneSolOneGrid;
	fSolHasNewGrid = FALSE;
	fIterOfGeneration = 0;
	fine = new TGrid( bt, this, input );
	if ( !fine ) NewtonFatalError( "memory allocation of TGrid failed" );
	grid = fine;
	coarse = new TGrid( bt, this, input );
	if ( !coarse ) NewtonFatalError( "memory allocation of TGrid failed" );
	fSolutionScaler = NewVector( variables );
	//	for ( int i = 0; i < variables; ++i ) {
	//		fSolutionScaler->vec[i] = 1.0;
	//	}
	auxmat = NewMatrix( variables, maxGridPoints, kColumnPointers );
    aux = NewVector( maxGridPoints + 2 );
    aux2 = NewVector( maxGridPoints + 2 );
    M = NewVector( maxGridPoints+1 );
	integral = New1DArray( variables );
	grid = fine;
	UpdateGridDimension();
	init = TRUE;	
//	check of values for grid correction
	if ( fEnd > fStart ) {
		if ( GetLeft() > fStart ) {
			cerr << "#warning: invalid value for 'start' in TAdaptiveGrid" << NEWL;
			fputc( '\a', stderr );
			fputc( '\a', stderr );
			fputc( '\a', stderr );
			fAlpha = 0.0;
		}
		if ( GetRight() < fEnd ) {
			cerr << "#warning: invalid value for 'end' in TAdaptiveGrid" << NEWL;
			fputc( '\a', stderr );
			fputc( '\a', stderr );
			fputc( '\a', stderr );
			fAlpha = 0.0;
		}
		if ( fAlpha >= 1.0 ) {
			cerr << "#warning: invalid value for 'end' in TAdaptiveGrid" << NEWL;
			fputc( '\a', stderr );
			fputc( '\a', stderr );
			fputc( '\a', stderr );
			fAlpha = 0.0;
		}
	}
}

TAdaptiveGrid::~TAdaptiveGrid( void )
{
	Free1DArray( integral );
	DisposeVector( M );
	DisposeVector( aux2 );
	DisposeVector( aux );
	DisposeMatrix( auxmat );
	DisposeVector( fSolutionScaler );
	delete coarse;
	delete fine;
}

void TAdaptiveGrid::switch_to_coarse( TBVPSolverPtr solver )
{
	int	gridPoints = coarse->GetNGridPoints();	
   
    grid = coarse;
	SetLeft( fine->GetLeft() );
	SetRight( fine->GetRight() );
	solver->UpdateAllDimensions( gridPoints );
}

void TAdaptiveGrid::switch_to_fine( TBVPSolverPtr solver )
{
	int	gridPoints = fine->GetNGridPoints();	
   
	grid = fine;
	solver->UpdateAllDimensions( gridPoints );
}

void TGrid::InitGrid( TNewtonPtr bt, TAdaptiveGridPtr /*adapGrid*/, TBVPInputPtr input )
{
	int					maxGridPoints = bt->GetMaxGridPoints();
	int 				variables = bt->GetNVariables();
	int 				initialGridPoints = bt->GetInitialGridpoints();
	
    x = NewVector( maxGridPoints );
    y = NewMatrix( variables, maxGridPoints, kColumnPointers );
	bcFlagLeft = new int[variables];
	if ( !bcFlagLeft ) NewtonFatalError( "memory allocation of TGrid failed" );
	bcFlagRight = new int[variables];
	if ( !bcFlagRight ) NewtonFatalError( "memory allocation of TGrid failed" );
    fBcLeft = NewVector( variables );
    fBcRight = NewVector( variables );
    yleft = NewVector( variables );
    yright = NewVector( variables );
	fFDWeightsConvM = NewVector( maxGridPoints );
	fFDWeightsConv = NewVector( maxGridPoints );
	fFDWeightsDiffM = NewVector( maxGridPoints );
	fFDWeightsDiff = NewVector( maxGridPoints );

	left	= input->fLeft;
	right	= input->fRight;

	int		*bcFLI = input->bcFlagLeft;
	int		*bcFRI = input->bcFlagRight;
	Double	*bcLI = input->fBcLeft->vec;
	Double	*bcRI = input->fBcRight->vec;
	Double	*bcL = fBcLeft->vec;
	Double	*bcR = fBcRight->vec;
	Double	*yLI = input->yleft->vec;
	Double	*yRI = input->yright->vec;
	Double	*yL = yleft->vec;
	Double	*yR = yright->vec;
	for ( int i = 0; i < variables; ++i ) {
		bcFlagLeft[i] = bcFLI[i];
		bcFlagRight[i] = bcFRI[i];
		bcL[i] = bcLI[i];
		bcR[i] = bcRI[i];
		yL[i] = yLI[i];
		yR[i] = yRI[i];
	}
	
	if ( initialGridPoints == 0 || initialGridPoints > maxGridPoints ){
    	n_gridpoints = maxGridPoints;
	}
	else{
		n_gridpoints = initialGridPoints;
	}

//	Make_equi_Grid( adapGrid );
}

TGrid::~TGrid( void )
{
	DisposeVector( fBcRight );
    DisposeVector( fBcLeft );
	delete bcFlagRight;
	delete bcFlagLeft;
   	DisposeMatrix( y );
    DisposeVector( x );
}

void TGrid::Make_equi_Grid( void )
{
    /*  generates equally spaced grid with N nodes, where
        'left' and 'right' are not part of the grid  */

    int         i;
    Double   	*locX = x->vec;
	Double		dx = ( right - left ) / ( n_gridpoints + 1 );

    if ( dx <= 0.0 ) {  /* errorchecking */
		NewtonFatalError( "grid generation failed: check boundaries" );	
	}

    for ( i = 1, locX[0] = left + dx; i < n_gridpoints; ++i ){
        locX[i] = locX[0] + i * dx;
    }
}

void TGrid::UpdateOneGridDimension( int number )
{
	x->len= number;
	y->cols= number;
}

void TGrid::ComputeFDWeights( void )
{
	Double	h, hm, hnenn;
	
	for ( int k = 1; k < n_gridpoints-1; ++k ) {
		h = x->vec[k+1] - x->vec[k];
		hm = x->vec[k] - x->vec[k-1];
		hnenn = h * hm * ( h + hm );
		fFDWeightsConvM->vec[k] = hm * hm  / hnenn;
		fFDWeightsConv->vec[k]  = h * h / hnenn;
		fFDWeightsDiffM->vec[k] = hm / hnenn;
		fFDWeightsDiff->vec[k]  = h / hnenn;
	}

}
