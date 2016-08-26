/*
	FLScan.h: Header file for the flamelet file scanner.
		
	© Josef Gšttgens, Peter Terhoeven, Ian Herwono, 1993
	
	Version 1.b2		07/21/93
*/

#ifndef __FLScan__
#define __FLScan__


#define YYOVERFLOW


#ifdef __cplusplus
extern "C" {
#endif

/*function prototypes (intern, used from flscan.main.c, flscan.y ) */
void InitFloatBuffers(int);
void ResetFloatBuffer(void);
void AddNumber(Double num);
void PrintFloatBuffer(void);
void PutSymH (char *id,char *unit,int type, int vali,Double valf,char *vals); 
void PutSymBody(char *id,char *unit,Double *vector,int len); 
#ifdef __cplusplus
}
#endif


/*
	Globals
*/
extern FILE *gFlex;					/*  diagnostic output for Flex			 */
extern FILE *gBison;				/*  diagnostic output for Bison 		 */
extern FILE *gFloatBuffer;			/*  diagnostic output for floatbuffer	 */


enum { kInitialState, kInHeader, kInBody, kInTrailer };

extern int	gSection;				/*	signals context					*/
extern int	gExpectArray;			/*	may be used to skip newlines	*/
extern int	gFLLine;					/*	keeps track of line numbers		*/

extern int 		gFBElems;
extern Double 	*gFBufhead;

#endif	/* __FLScan__ */
