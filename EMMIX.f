      PROGRAM EMMIX 

C      Version 1.3 1999 

C   PURPOSE 

C      The main purpose of this program is to fit a mixture model of
C      multivariate normal or t-distributions to a user supplied data set
C      This is done via the EM algorithm.  A large number of other features
C      are included that were found to be of use when fitting mixture models.
C
C      For information about how to run and use this program see 
C      http:\\www.maths.uq.edu.au\~gjm\emmix\emmix.html 

C   HISTORY

C      D.Peel Nov 1995 (Called NMM)
C      Combined with MMresamp (D.Peel Nov 1995) on Oct 1996
C      Renamed MIXCLUS in 1996
C      Renamed MIXFIT May 1997
C      Renamed EMMIX from MIXFIT version 1.3 Oct 1998
C      Version EMMIX 1.2     
C      Version EMMIX 1.3 fixed bugs from previous version and modified
C                        structure to be used as a DLL in other versions
C                        May 1999

C     This is the header off the program and calls
C     the main loop of the program + the user interface

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX.max'

C    Global Parameters
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO

C   INPUT PARAMETERS
      INTEGER   FLAGS(40),     ! Flags summarising options chosen (see below)
     &          FYLENO,        ! File id for output file
     &          NIND,          ! Number of points (samples, or observations)
     &          NATT,          ! Number of dimensions (variates or attributes)
     &          NG0,           ! Minimum number of groups to be tested
     &          NG1,           ! Maximum number of groups to be tested
     &          NCOV,          ! Structure of covariance matrices
     &          IX,IY,IZ       ! Random seeds

      DIMENSION X(MNIND,MNATT) ! Data or sample. 

C   OTHER INPUT PARAMETERS
      INTEGER   TIDT(MNIND),              ! User defined partition of sample
     &          USA(MNIND),               ! Grouping of classified sample
     &          RNATT,                    ! Reduced number of dimensions
     &          RV(MNATT)                 ! Array of dimensions to be used
      DIMENSION TXUU(MAXNG),              ! User defined NU for fitting t-dist
     &          TXVAR(MAXNG,MNATT,MNATT), ! User defined covariance matrices
     &          TT(MAXNG),                ! User defined mixing proportions
     &          TXMU(MAXNG,MNATT),        ! User defined group means
     &          W(MNIND,MAXNG),           ! User defined posterior probabilities
     &           TOLS(4)                  ! User stopping tolerances for EM

C   OUTPUT PARAMETERS
      INTEGER    IDTS(MNIND,MAXNG)                !Partitions of the data
      DIMENSION  XMUS(MAXNG,MAXNG,MNATT),         !Group means
     &           XVARS(MAXNG,MAXNG,MNATT,MNATT),  !Group covariance matrices
     &           XUU(MAXNG),                      !Group NU values (for t-dist)
     &           TS(MAXNG,MAXNG),                 !Group mixing proportions
     &           AIC(MAXNG),                      !Akaike Information Criterion
     &           BIC(MAXNG),                      !Bayesian  "  "       "  "
     &           AWE(MAXNG),                      !Approx. Weight. Evidence
     &           PVAL(MAXNG),                     !P-value of -2log(Lambda)
     &           TLL(MAXNG)                       !-2log(Lambda)
C   OTHER OUTPUT PARAMETERS
      DIMENSION  SET(MAXNG),            !Standard Errors of mixing proportions
     &           SEU(MAXNG,MNATT),      !Standard Errors of group means
     &           SEV(MAXNG,MNATT,MNATT) !Standard Errors of covariance matrices
 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   FLAGS CODES
   
C    1 - % of data used to form random starts (100 =std random start)
C    2 - SEM FLAG  (0-normal EM, 1-Stochastic EM)
C    3 - temp 1- tru data fit 2- bootstrap fit (no output to screen) 
C                                   3 -Bootstrap under H0
C    4 - Type of start 1-partition, 2-parameter 3-auto 4-weights
C    5 - Number of k-means starts 
C    6 - Display density values to use as a discriminant rule -disc
C    7 - T density (1-T ,0-normal)
C    8 - 0-simulate 1-Bootstrap analysis, 2-Specific analysis,
C        3-Full auto analysis, 4-Discriminant, 5- Prediction 
C    9 - 1-Final EM iterations / 2-Initial EM iterations
C   10 - resamp test (0-No, >0 -yes (Number of replications))
C   11 - Space efficient version (0-no 1- partial 2- minimal)
C   12 - Partial user allocation knowledge (0-no,1-yes)
C   13 - unused
C   14 - Weighted data set (0-no 1-yes) 
C   15 - Output Index+partition for external plot (0-no, 1=yes)
C   16 - Output boot distrib for external plot (0-no,1-yes)
C   17 - Estimate Standard Errors (0-no >0 = No its or =1 yes)
C   18 - SE's Method (0-parametric,1-Samp w/replace,2-weight like,4-info method 
C   19 - Variable Selection:-1 samp 1- adjust data 2- adjust parameters as well 
C   20 - Output to separate file  1- parameters,  2-point likelihoods 3-data
C   21 - Use Aitken's acceleration during bootstrapping (<0 active >0 on)
C   22 - Output subset of data to separate file 
C   23 - Read Parameters
C   24 - Read Partition
C   25 - Read Posterior
C   26 - Number of random starts


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   ERROR CODES

C   1 - Covariance matrix pivot zero (ie close to singular)
C   2 - Covariance matrix is not positive semi-definite
C   4 - Nullity = 0
C   5 - Determinant = 0
C   6 - Input partition incorrect
C   11 - Number of data points too big for this compilation
C   12 - Number of data variables too big for this compilation
C   13 - Unused 
C   14 - Maximum Number of clusters too big for this compilation
C   15 - Number of clusters too big for this compilation
C   21 - Not enough points in cluster at initial estimation stage 
C   22 - No points allocated to cluster during an EM iteration 
C   23 - Problem in the generation of a bootstrap sample
C   25 - Estimated Nu value when fitting T's is < or equal to Zero
C   31 - No stable starting solution could be found
C   40 - Random number generator not working 
C  -41 - Warning : k-means reached maximum number of iterations 
C  -53 - Warning : Estimated Nu value when fitting T's limited to 300 
C -111 - Warning : Some points have zero likelihood


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  INPUT/OUTPUT FILE ID NUMBERS 

C 21 - main data file + starting parameters or partition
C 22 - main output file from main gives clusterings
C 56 - optional allocation for export to external plotting package
C 57 - optional bootstrap   "   "   "   "       "         "    "
C 28 - 'hier.inp' optional input file specifies hierarchical methods
C 42 - 'respSE.out'
C 42 - 'respH0.out' output file for fit under H0 for last bootstrap replicate
C 42 - 'respH1.out' output file for fit under H1 for last bootstrap replicate
C 43 -  output file of bootstrap sample for last -2logLambda replicate
C 43 -  output file of bootstrap sample for last SE replicate
C 25 - 'boot_vs_.out output file contain bootstrap replicates of -2logL 
C 26 - Standard error estimates of parameters for replications 
C 29 - Output parameters for variable selection
C 39 - Output data when a subset of variables are used

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C    Check to see if the random number generator works (try two formats)
      CALL DETERRANDOM(IER)
 
C -------------------- Snip Here ----------------------------------
C    For DLL remove this section

C    Read in parameters and options from the user
      CALL SETUP(NIND,NATT,NG0,NG1,NCOV,X,TOLS,USA,
     &   RV,RNATT,SIG,TXUU,TXMU,TXVAR,TT,TIDT,W,IER)
      IF (IER.GT.0) THEN
c      ERROR as the input file may be in the wrong format of the 
c      specified parameters do not match the input file.
       GOTO 99
      ENDIF

C -------------------- End Snip -----------------------------------

C    Call the main control section of the program
      CALL MAIN(NIND,NATT,NG0,NG1,NCOV,X,TOLS,USA,
     &   RV,RNATT,SIG,TXUU,TXMU,TXVAR,TT,TIDT,W,
     &   XMUS,XVARS,TS,XUU,AIC,BIC,AWE,TLL,PVAL,IDTS,SET,SEU,SEV,IER)
      
      WRITE (*,*) 'Program run complete'
99    CONTINUE
      END 
C
C
C
      SUBROUTINE SETUP(NIND,NATT,NG0,NG1,NCOV,X,TOLS,
     &         USA,RV,RNATT,SIG,XUU,TXMU,TXVAR,TT,TIDT,W,IER)
C      Read in main parameters and set options with interactive
C      questions with the user (or alternatively from a file)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,THING,RNATT,RV(MNATT),POINT,FI
      INTEGER USA(MNIND)
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO
      COMMON /STORE3/ OFYLE,INFYLE,P1FYLE,P2FYLE,P3FYLE,OFYLE2
      DIMENSION X(MNIND,MNATT),TOLS(4),
     &          TX(MNIND,MNATT),XUU(MAXNG)
      DIMENSION TXMU(MAXNG,MNATT,MNATT),TXVAR(MAXNG,MNATT,MNATT),
     &          TT(MAXNG),W(MNIND,MAXNG)
      INTEGER   TIDT(MNIND)
      CHARACTER INFYLE*20,OFYLE*20,ANSWER*20,OFYLE2*20
      CHARACTER P1FYLE*20,P2FYLE*20,P3FYLE*20

C     Set all FLAGS to default zero
      DO 1 I=1,40
       FLAGS(I)=0
1     CONTINUE

      TOLS(1)=TAUTO
      TOLS(2)=MITAUT
      TOLS(3)=TFINAL
      TOLS(4)=MITFIN

      WRITE (*,501)
501    FORMAT (///////////,
     &/,12X,'------------------------------------------------------',
     &/,15X,'       _____  __   __    __   __ __ _    _ ',
     &/,15X,'      | ____| | \\_/ |    | \\_/ | || \\\\  // ',
     &/,15X,'      ||____  ||\\_/||    ||\\_/|| ||  \\\\//  ',
     &/,15X,'      | ____| ||   || __ ||   || ||   ||  ',
     &/,15X,'      ||____  ||   || -- ||   || ||  //\\\\  ',
     &/,15X,'      |_____| ||   ||    ||   || || //  \\\\ ',
     &/,15X,' ' ,
     &/,12X,'------------------------------------------------------',
     &/,15X,'          EM based MIXTURE program' ,
     &/,15X,'               Version 1.3  1999              ',
     &/,12X,'------------------------------------------------------')
503   CONTINUE
      WRITE(*,519)
519   FORMAT(
     &/,12X,' Do you wish to:',
     &/,12X,'  0. Simulate a sample from a normal mixture model',
     &/,12X,'  1. Carry out a bootstrap-based assessment of',
     &/,12X,'     standard errors and/or the number of components (g)',
     &/,12X,'  2. Fit a g-component normal mixture model for a',
     &/,12X,'     specified g',
     &/,12X,'  3. Fit a g-component normal mixture model for a',
     &/,12X,'     range of values of g',
     &/,12X,'  4. Perform discriminant analysis',
     &/,12X,'  5. Make predictions for new data',
     &/,12X,'  6. Form parameter estimates from data + allocation',
     &/,12X,'------------------------------------------------------',
     &/)
      READ (*,*)  FLAGS(8)
      IF ((RANDTYPE.EQ.0).AND.(FLAGS(8).EQ.1)) THEN
        WRITE (*,*) 'PROBLEM: Unable to use Bootstrap/Standard Error'
        WRITE (*,*) '         since the random number does not seem'
        WRITE (*,*) '         to work.  Only options which do not use'
        WRITE (*,*) '         random numbers are available'
        GOTO 503
      ENDIF

      IF (FLAGS(8).LE.1) FLAGS(23)=1
      IF (FLAGS(8).EQ.5) FLAGS(23)=1
      IF (FLAGS(8).EQ.6) FLAGS(24)=1


      WRITE (*,*)'Enter name of input file: '
      READ (*,'(A)') INFYLE
      OPEN (UNIT=21,FILE=INFYLE,STATUS = 'OLD',ERR=505)
      GOTO 506

505   WRITE (*,*) 'Cannot locate ',INFYLE
      WRITE (*,*) 'please re-enter file name: '
      READ (*,'(A)') INFYLE
      OPEN (UNIT=21,FILE=INFYLE,STATUS = 'OLD',ERR=505)

506   CONTINUE
      IF (FLAGS(8).NE.1) THEN
        WRITE (*,*)'Enter name of output file: '
        READ (*,'(A)') OFYLE
        OPEN (UNIT=22,FILE=OFYLE,STATUS = 'UNKNOWN')
      ENDIF
      THING=5
      WRITE(*,*)
      IF (FLAGS(8).EQ.1) THEN
        WRITE(*,*)'Do you want:'
        WRITE(*,*)'  1. A Bootstrap analysis of -2log(Lambda)'
        WRITE(*,*)'  2. A Standard Error analysis'
        WRITE(*,*)'  3. Both 1 and 2'
        READ (*,*) METHOD
        IF ((METHOD.EQ.1).OR.(METHOD.EQ.3)) THEN
          WRITE (*,*)' Enter name of output file for Bootstrap: '
          READ (*,'(A)') OFYLE
          OPEN(UNIT=25,FILE=OFYLE,STATUS='UNKNOWN')
        ENDIF
        IF ((METHOD.EQ.2).OR.(METHOD.EQ.3)) THEN
          WRITE (*,*)' Enter name of output file for',
     &    ' Standard Errors: '
          READ (*,'(A)') OFYLE
          OPEN(UNIT=26,FILE=OFYLE,STATUS='UNKNOWN')
        WRITE (*,*)'Which method of estimation: '
        WRITE (*,*)' 0 Parametric'
        WRITE (*,*)' 1 Sampling with replacement'
        WRITE (*,*)' 2 weighted likelihood'
        WRITE (*,*)' 4 information based method '
        WRITE (*,*)'   (Only available unequal covariances)'
          READ (THING,*) FLAGS(18)
          IF (FLAGS(18).GE.1) THEN
            IF (METHOD.EQ.3) THEN
      WRITE(*,*)'Please place the original sample at the start of'
      WRITE(*,*)'the input file'
            ELSE
      WRITE(*,*) 'Please put original sample to the input file '
            ENDIF
          ENDIF
        ENDIF
        IF (METHOD.EQ.1) THEN
          NEDRAN=1
          WRITE(*,*)' [Warning may need extensive time]'
          WRITE(*,*)'  How many bootstrap replications'
          READ (*,*) FLAGS(10)
          FLAGS(17)=0
        ELSEIF ((METHOD.EQ.2).AND.(FLAGS(18).LT.4)) THEN
          NEDRAN=1
          WRITE(*,*)' [Warning may need extensive time]'
          WRITE(*,*)'  How many replications to estimate'
          WRITE(*,*)'  the Standard Errors'
          READ (*,*) FLAGS(17)
          FLAGS(10)=0
        ELSEIF (METHOD.EQ.3) THEN
          NEDRAN=1
          WRITE(*,*)' [Warning may need extensive time]'
          WRITE(*,*)'  How many replications to estimate'
          WRITE(*,*)'  the Standard Errors/Bootstrap'
          READ (*,*) FLAGS(10)
          FLAGS(17)=FLAGS(10)
        ELSEIF (FLAGS(18).EQ.4) THEN
          FLAGS(17)=1
        ENDIF

      ELSEIF (FLAGS(8).EQ.3) THEN
        WRITE(*,*)'Do you wish to carry out a bootstrap test'
        WRITE(*,*)' to assess the number of components (Yes/No)- '
        READ (THING,'(A)') ANSWER
        IF   ((ANSWER.EQ.'yes').OR.(ANSWER.EQ.'YES')
     &    .OR.(ANSWER.EQ.'y').OR.(ANSWER.EQ.'Y')) THEN
          IF (RANDTYPE.EQ.0) THEN
            WRITE (*,*) 'Unable to use bootstrap as random number'
            WRITE (*,*) ' generator does not work'
            FLAGS(10)=0
          ELSE
            NEDRAN=1
            WRITE(*,*)' [Warning may need extensive time]'
            WRITE(*,*)'  How many bootstrap replications: '
            READ (*,*) FLAGS(10)
          ENDIF
        ELSE
          FLAGS(10)=0
        ENDIF
      ENDIF

      WRITE(*,*) 'Number of entities: '
      READ (THING,*) NIND
      IF (NIND.GT.MNIND) THEN
       WRITE (*,*)
       WRITE (*,*) 'ERROR: number of entities too large modify MNIND '
       WRITE (*,*) 'parameter to ',NIND,' in file EMMIX.max'
       WRITE (*,*) 'and recompile.'
       IER=11
       GOTO 599
      ENDIF
      IF ((NIND.GT.1000).AND.(FLAGS(8).NE.0)) THEN
       WRITE (*,*) '[Maybe you should consider selecting the space'
       WRITE (*,*) 'saving version from the other options]'
      ENDIF
      IF ((FLAGS(8).NE.0).OR.(FLAGS(8).NE.1)) THEN
       WRITE(*,*) 'Total Number of variables/dimensions'
       WRITE(*,*) ' in the input file: '
       READ (THING,*) NATT
       WRITE (*,1505) NATT
1505   FORMAT(' How many variables to be used in the analysis',/,
     &   ' (re-enter',I3,' if you wish to use all the variables): ')
c       WRITE (*,*)'How many variables in the subset to be used'
c       WRITE (*,*)' (re-enter ',NATT,' if you wish to use all'
c       WRITE (*,*)'   the variables): '
       READ (*,*) RNATT
       IF (RNATT.NE.NATT) THEN
        WRITE (*,*)'Please type the variables that are to be '
        WRITE (*,*)'included (space delimited on a single line): '
        READ (*,*) (RV(II),II=1,RNATT)
        IF (FLAGS(8).LE.1) THEN
         FLAGS(19)=-1
        ELSE
         FLAGS(19)=1
         WRITE (*,*) 'Do you wish to have the data subset written'
         WRITE (*,*) 'written to file (0-no, 1-yes): '
         READ (*,*) TEMP
         IF (TEMP.EQ.1) THEN
          FLAGS(22)=1
          WRITE (*,*) 'What do you wish this file to be called: '
          READ (THING,'(A)') ANSWER
          OPEN(UNIT=39,FILE=ANSWER,STATUS='UNKNOWN')
         ENDIF
       ENDIF
       IF (NATT.GT.MNATT) THEN
       WRITE(*,*)'ERROR: number of variables too large modify MNATT'
       WRITE (*,*) 'parameter to ',RNATT,' in file EMMIX.max '
       WRITE (*,*) 'and recompile.'
       IER=12
       GOTO 599
       ENDIF
       ELSE
      IF (NATT.GT.MNATT) THEN
       WRITE (*,*)
       WRITE(*,*)'ERROR: number of variables too large modify MNATT'
       WRITE (*,*) 'parameter to ',NATT,' in file EMMIX.max '
       WRITE (*,*) 'and recompile.'
       IER=12
       GOTO 599
      ENDIF
      ENDIF
      ELSE
       WRITE(*,*) 'Number of variables/dimensions: '
       READ (THING,*) NATT
      IF (NATT.GT.MNATT) THEN
       WRITE (*,*)
       WRITE(*,*)'ERROR: number of variables too large modify MNATT'
       WRITE (*,*) 'parameter to ',NATT,' in file EMMIX.max '
       WRITE (*,*) 'and recompile.'
       IER=12
       GOTO 599
      ENDIF
      ENDIF

      IF (FLAGS(8).EQ.3) THEN
502    WRITE(*,*)'What is the minimum number of components',
     & ' you wish to test (eg 1): '
       READ (THING,*) NG0
       WRITE(*,*) 'What is the maximum number of components',
     & ' you wish to test (eg 10): '
       READ (THING,*) NG1
       IF (NG0.GT.NG1) THEN
         WRITE (*,*) '   Min must be less than max try again'
         GOTO 502
       ENDIF
       IF (FLAGS(10).GT.0) THEN
        WRITE (*,*)' Do you wish to stop when P-value is',
     &   ' insignificant (0-No,1-Yes): '
        READ (*,*) TEMP
        IF (TEMP.EQ.1) THEN
         WRITE (*,*) 'What level of significance (ie 10 =10%): '
         READ (*,*) SIG
        ENDIF
       ENDIF
      ELSEIF (FLAGS(8).EQ.0) THEN
       NEDRAN=1
       WRITE(*,*) 'How many components do you want to generate: '
       READ (THING,*) NG0
       NG1=NG0

      ELSEIF ((FLAGS(8).EQ.2).OR.(FLAGS(8).EQ.4)
     &        .OR.(FLAGS(8).EQ.5).OR.(FLAGS(8).EQ.6)) THEN
       WRITE(*,*) 'How many components do you want to fit: '
       READ (THING,*) NG0
       NG1=NG0

      ELSEIF (FLAGS(10).GT.0) THEN
      WRITE(*,*)'What value of g do you wish to test (g vs g+1): '
       READ (THING,*) NG0
       NG1=NG0+1
      ELSEIF (FLAGS(17).GT.0) THEN
       WRITE(*,*) 'How many components are you fitting: '
       READ (THING,*) NG0
       NG1=NG0
      ENDIF
      IF (NG0.GT.MAXNG) THEN
        WRITE (*,*)
        WRITE(*,*)'ERROR: number of components too large modify MAXNG'
        WRITE (*,*) 'parameter to ',NG0,' in file EMMIX.max '
        WRITE (*,*) 'and recompile.'
        IER=15
        GOTO 599
      ELSEIF (NG1.GT.MAXNG) THEN
        WRITE (*,*)
        WRITE(*,*)'ERROR: upper number of components too large modify'
        WRITE (*,*) 'MAXNG parameter to ',NG1,' in file EMMIX.max'
        WRITE (*,*) 'and recompile.'
        IER=14
        GOTO 599
      ENDIF

      IF (FLAGS(8).EQ.0) THEN
       GOTO 1010
      ENDIF

      IF (NG1.GT.1) THEN
510   WRITE(*,*)'Covariance matrix option (1 = equal,2 = unrestricted,'
      WRITE(*,*)'   3 = diagonal equal,4 = diagonal unrestricted): '
       READ (THING,*) NCOV
       IF ((THING.EQ.5).AND.((NCOV.LT.1).AND.(NCOV.GT.4))) THEN
         WRITE (*,*) '  ERROR expecting a 1,2,3 or 4 repeat answer: '
         GOTO 510
       ENDIF

       IF (FLAGS(8).EQ.6) GOTO 597

       IF (FLAGS(8).EQ.2) THEN
520      WRITE(*,*) 'Switch for initialisation'
         WRITE(*,*)'   (1 = initial outright grouping,'
         WRITE(*,*)'    2 = initial parameter estimates,'
         WRITE(*,*)'    3 = automatic initial grouping '
         WRITE(*,*)'    4 = initial soft or fractional grouping): '
         READ (THING,*) FLAGS(4)
        IF ((THING.EQ.5).AND.(FLAGS(4).NE.1).AND.
     &                  (FLAGS(4).NE.2).AND.(FLAGS(4).NE.4)
     &                  .AND.(FLAGS(4).NE.3)) THEN
         WRITE (*,*) '  ERROR expecting a 1,2 or a 3 repeat answer: '
         GOTO 520
        ENDIF

        IF (FLAGS(4).EQ.1) FLAGS(24)=1
        IF (FLAGS(4).EQ.2) FLAGS(23)=1
        IF (FLAGS(4).EQ.4) FLAGS(25)=1

         IF ((FLAGS(19).EQ.1).AND.((FLAGS(4).EQ.2).OR.
     &      (FLAGS(8).EQ.5))) THEN
          WRITE (*,*) ' Are the parameters for the '
          WRITE(*,*) '  1. variables subset'
          WRITE(*,*) '  2. original variables: '
          READ (*,*) FLAGS(19)
         ENDIF
       ELSEIF (FLAGS(8).EQ.5) THEN
         FLAGS(4)=2
         FLAGS(6)=1
       ELSEIF (FLAGS(8).EQ.1) THEN
         FLAGS(4)=3
       ELSEIF (FLAGS(8).EQ.4) THEN
           FLAGS(4)=1
           FLAGS(6)=2
           FLAGS(12)=1
       ELSE
         FLAGS(4)=3
       ENDIF
       IF ((FLAGS(4).EQ.3).AND.(METHOD.NE.2)) THEN
         IF (HIRFLG.EQ.0) THEN
          WRITE(*,*) '  This version has no hierarchical methods'
          WRITE(*,*) '  if needed recompile with HIRFLG=1'
         ELSEIF (NIND.GT.400) THEN
          WRITE(*,*) '  This is a lot of data entities to run a'
          WRITE(*,*) '  hierarchical method on it may be time'
          WRITE(*,*) '  consuming -see readme file'
         ENDIF
         WRITE(*,*)'How many random starts: '
         READ (THING,*) FLAGS(26) 
          IF ((RANDTYPE.EQ.0).AND.(FLAGS(26).GT.0)) THEN
            WRITE (*,*) 'Unable to use random starts as random'
            WRITE (*,*) 'number generator does not work'
            FLAGS(26)=0
          ELSE
            NEDRAN=1
          ENDIF
          IF (FLAGS(26).GT.0) THEN
c           WRITE(*,*)'Do you want to use'
c           WRITE(*,*)' 1- Standard random starts'
c           WRITE(*,*)' 2- Subset random starts: '
c           READ (THING,*) FLAGS(1)
c           IF (FLAGS(1).EQ.1) THEN
c             FLAGS(1)=0
c           ELSE
1009        WRITE(*,*)'What percentage of the data is to'
            WRITE(*,*)'be used to form random starts: '
            READ (THING,*) FLAGS(1)
            IF ((FLAGS(1).GT.100).OR.(FLAGS(1).LT.0)) THEN
             WRITE(*,*) 'Error: This answer must be a percentage'
             WRITE(*,*) '       (ie between 0 and 100)'
             GOTO 1009
            ELSEIF (FLOAT(FLAGS(1))/100*FLOAT(NIND).LE.(NATT+1)) THEN
             WRITE(*,*) 'Warning: There may be not enough points'
             WRITE(*,*) '         used in the random starts resulting'
             WRITE(*,*) '         in singular covariance matrices'
            ENDIF
          ENDIF
         WRITE(*,*)'How many k-means starts: '
         READ (THING,*) FLAGS(5)
          IF ((RANDTYPE.EQ.0).AND.(FLAGS(5).GT.0)) THEN
            WRITE (*,*) 'Unable to use k-means as random'
            WRITE (*,*) 'number generator does not work'
            FLAGS(5)=0
          ELSE
            NEDRAN=1
          ENDIF
       ENDIF
      ENDIF
       IF (FLAGS(8).EQ.1) THEN
         FLAGS(4)=3
       ENDIF
CCCCCCCc      IF (FLAGS(8).LE.4) THEN
      IF (FLAGS(8).LE.5) THEN
       WRITE(*,*)'Are extra options required(Y/N): '
       READ (THING,'(A)') ANSWER
       IF   ((ANSWER.EQ.'yes').OR.(ANSWER.EQ.'YES')
     &  .OR.(ANSWER.EQ.'y').OR.(ANSWER.EQ.'Y')) THEN
      CALL EXOPT(NG1,TOLS,XUU,NCOV,NEDRAN)
      ENDIF
      ENDIF

1010    IF (NEDRAN.EQ.1) THEN
        IF (RANDTYPE.EQ.0) THEN
          WRITE (*,*) 'WARNING: You have selected options which'
          WRITE (*,*) '         utilise random numbers but the'
          WRITE (*,*) '         random number generator is not'
          WRITE (*,*) '         working'
        ELSEIF (RANDTYPE.EQ.1) THEN
          WRITE(*,*)'Random seeds 3 seeds needed : '
          WRITE(*,*)'  random seed 1 [0-30000]: '
          READ (THING,*) IX
          WRITE(*,*)'  random seed 2 [0-30000]: '
          READ (THING,*) IY
          WRITE(*,*)'  random seed 3 [0-30000]: '
          READ (THING,*) IZ
        ENDIF
      ENDIF
      WRITE(*,535)
535    FORMAT (////////////////////////)
        FI=6
        CALL SUMRY(NIND,NATT,NG0,NG1,NCOV,TOLS,NEDRAN,FI,
     &             RV,RNATT)
        IF (FLAGS(8).EQ.1) THEN
          IF (FLAGS(17).GT.0) THEN
             FI=26
        CALL SUMRY(NIND,NATT,NG0,NG1,NCOV,TOLS,NEDRAN,FI,
     &             RV,RNATT)
          ENDIF
          IF (FLAGS(10).GT.0) THEN
            FI=25
        CALL SUMRY(NIND,NATT,NG0,NG1,NCOV,TOLS,NEDRAN,FI,
     &             RV,RNATT)
          ENDIF
        ELSE
          FI=22
        CALL SUMRY(NIND,NATT,NG0,NG1,NCOV,TOLS,NEDRAN,FI,
     &             RV,RNATT)
        ENDIF

C      Output to screen the form of the input file required
        CALL INPSUM(NIND,NATT,NG0)

597     CONTINUE


      IF ((FLAGS(8).GT.1).OR.(FLAGS(18).GE.1)) THEN
cc     Modification D.Podlich 1995
        READ (21,*) ((X(I,J),J=1,NATT),I=1,NIND)
        WRITE (*,*) ' Read in data.'
      ENDIF

cc        Old form of data read
cc       DO 598 I=1,NIND
cc         READ (21,*) (X(I,J),J=1,NATT)
cc598    CONTINUE

      IF (FLAGS(19).GT.0) THEN
        DO 710 I=1,NIND
          DO 700 J=1,RNATT
           TX(I,J)=X(I,RV(J))
700       CONTINUE
          DO 710 J=1,RNATT
           X(I,J)=TX(I,J)
710     CONTINUE
        IF (FLAGS(19).EQ.1) THEN
          NATT=RNATT
        ENDIF
      ENDIF

      DO 111 II=1,NIND
       USA(II)=-1
111   CONTINUE

C     Setup partial user allocation
      IF (FLAGS(12).GT.0) THEN
16     CONTINUE
        READ (21,*) POINT,GRP
        IF (POINT.EQ.(-1)) GOTO 17
        USA(POINT)=GRP
       GOTO 16
17    CONTINUE
      ENDIF

       IF (FLAGS(8).EQ.1) THEN
          IF (FLAGS(17).GT.0) THEN
             FYLENO=26
          ENDIF
          IF (FLAGS(10).GT.0) THEN
            FYLENO=25
          ENDIF
       ELSE
          FYLENO=22
       ENDIF

      IF (FLAGS(23).EQ.1) THEN
       CALL USRPARAMETERS(NIND,NATT,NG0,NCOV,TXMU,TXVAR,
     &           TT,RNATT,RV)
      ENDIF

       IF (FLAGS(24).EQ.1) THEN
        CALL USRALLOC(NIND,NATT,NG0,TIDT,IER)
        IF (IER.EQ.6) RETURN
       ENDIF

       IF (FLAGS(25).EQ.1) THEN
C       Read in posterior probabilities from input file
        DO 24 I=1,NIND
         READ(21,*) (W(I,J),J=1,NG0)
24      CONTINUE
         write (*,*) ' Read in Posterior Probs.'
C       Correct posterior probabilities for classified data
        CALL PARCORR(NIND,NATT,NG0,USA,W)
       ENDIF

599    RETURN
      END



       SUBROUTINE EXOPT(NG1,TOLS,XUU,NCOV,NEDRAN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,THING,TEMP
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO
      COMMON /STORE3/ OFYLE,INFYLE,P1FYLE,P2FYLE,P3FYLE,OFYLE2
      DIMENSION TOLS(4),XUU(MAXNG)
      CHARACTER INFYLE*20, OFYLE*20, OFYLE2*20
      CHARACTER P1FYLE*20,P2FYLE*20,P3FYLE*20
        THING=5
700     CONTINUE
        WRITE (*,*) '              EXTRA OPTIONS       '
        WRITE(*,*)' ---------------------------------------------'
        WRITE (*,*) 'Please select option (selection will toggle):'
        IF (FLAGS(2).EQ.1) THEN
         WRITE(*,*)'  1. Stochastic EM option : YES'
        ELSE
         WRITE(*,*)'  1. Stochastic EM option : NO'
        ENDIF
        WRITE(*,*)'  2. Modify EM stopping criteria'
        IF (FLAGS(11).EQ.2) THEN
        WRITE(*,*)'  3. Space efficiency : EXTREME'
        ELSEIF (FLAGS(11).EQ.1) THEN
        WRITE(*,*)'  3. Space efficiency : MODERATE'
        ELSE
        WRITE(*,*)'  3. Space efficiency : OFF'
        ENDIF
        IF ((FLAGS(15).GT.0).OR.(FLAGS(16).GT.0)
     &   .OR.(FLAGS(20).GT.0)) THEN
        WRITE(*,*)'  4. Change extra output files'
        ELSEIF (FLAGS(8).EQ.1) THEN
        WRITE(*,*)'  4. Add extra output files : N/A'
        ELSE
        WRITE(*,*)'  4. Add extra output files'
        ENDIF
cccccc        IF (FLAGS(8).EQ.1) THEN
        IF (FLAGS(8).LE.1) THEN
        WRITE(*,*)'  5. Partial user classified data : N/A'
        ELSEIF (FLAGS(12).EQ.0) THEN
        WRITE(*,*)'  5. Partial user classified data : NO'
        ELSE
        WRITE(*,*)'  5. Partial user classified data : YES'
        ENDIF
        IF (FLAGS(7).GT.0) THEN
        WRITE(*,*)'  6. Estimate Standard Errors : N/A'
        FLAGS(17)=0
        ELSEIF (FLAGS(17).GT.1) THEN
        WRITE(*,*)'  6. Estimate Standard Errors : YES'
        ELSEIF (FLAGS(17).EQ.0) THEN
        WRITE(*,*)'  6. Estimate Standard Errors : NO'
        ENDIF
        WRITE(*,*)'  7. Fit mixture of Factor analysers : N/A '
      IF (((FLAGS(4).EQ.2).OR.(NG1.EQ.1)).AND.(FLAGS(8).NE.1)) THEN
        IF (FLAGS(6).EQ.0) THEN
        WRITE(*,*)'  8. Display discriminant density values : NO'
        ELSE
        WRITE(*,*)'  8. Display discriminant density values : YES'
        ENDIF
        ELSE
        WRITE(*,*)'  8. Display discriminant density values : N/A'
        ENDIF
        IF (FLAGS(10).GT.0) THEN
        WRITE(*,*)'  9. Change component distributions : N/A'
        ELSE
        WRITE(*,*)'  9. Change component distributions'
        IF (FLAGS(7).EQ.0) THEN
        WRITE(*,*)'     (Currently fitting NORMAL components)'
        ELSEIF (FLAGS(7).EQ.1) THEN
        WRITE(*,*)'     (Currently fitting t components'
        WRITE(*,*)'     with user defined NUs)'
        ELSEIF (FLAGS(7).EQ.2) THEN
        WRITE(*,*)'     (Currently fitting t components with'
        WRITE(*,*)'     iteratively calculated NU started at user'
        WRITE(*,*)'     estimate with unequal NUs for all groups)'
        ELSEIF (FLAGS(7).EQ.3) THEN
        WRITE(*,*)'     (Currently fitting t components with'
        WRITE(*,*)'     iteratively calculated NU started at user'
        WRITE(*,*)'     estimates with equal NUs for all groups'
        ELSEIF (FLAGS(7).EQ.4) THEN
        WRITE(*,*)'     (Currently fitting t components with'
        WRITE(*,*)'     iteratively calculated NU started at moment'
        WRITE(*,*)'     estimates with unequal NUs for all groups'
        ENDIF
        ENDIF
        IF ((FLAGS(10).GT.0).AND.(FLAGS(21).NE.1)) THEN
        WRITE(*,*)' 10. Use Aitken acceleration when bootstrapping'
        WRITE(*,*)'     -2log(lambda): NO'
        ELSEIF ((FLAGS(10).GT.0).AND.(FLAGS(21).EQ.1)) THEN
        WRITE(*,*)' 10. Use Aitken acceleration when bootstrapping'
        WRITE(*,*)'     -2log(lambda): YES'
        ELSE
        WRITE(*,*)' 10. Use Aitken acceleration when bootstrapping'
        WRITE(*,*)'     -2log(lambda): N/A '
        ENDIF
        WRITE(*,*)'  0. Run program '
        WRITE(*,*)' ---------------------------------------------'
        WRITE(*,*)
        READ (THING,*) TEMP
        IF (TEMP.EQ.1) THEN
          FLAGS(2)=1-FLAGS(2)
         IF ((RANDTYPE.EQ.0).AND.(FLAGS(2).NE.0)) THEN
          WRITE (*,*) 'Unable to use Stochastic EM as random'
          WRITE (*,*) 'number generator does not work'
          FLAGS(2)=0
         ELSE
          NEDRAN=1
         ENDIF
        ELSEIF (TEMP.EQ.2) THEN
       WRITE(*,*)'  -Set tolerance automatic methods (Default=',
     &               TAUTO,')'
       WRITE(*,*)'   Either set new value or 0 for default: '
          READ (THING,*) TOLS(1)
          IF (TOLS(1).EQ.0) THEN
           TOLS(1)=TAUTO
          ENDIF
        WRITE(*,*)'  -Set max number of iterations for automatic'
        WRITE(*,*)'   methods (Default=',MITAUT,')'
        WRITE(*,*)'    Either set new value or 0 for default: '
          READ (THING,*) TOLS(2)
          IF (TOLS(2).EQ.0) THEN
           TOLS(2)=MITAUT
          ENDIF
          WRITE(*,*)'  -Set tolerance final fit (Default=',
     &                      TFINAL,')'
          WRITE(*,*)'   Either set new value or 0 for default: '
          READ (THING,*)  TOLS(3)
           IF (TOLS(3).EQ.0) THEN
             TOLS(3)=TFINAL
           ENDIF
       WRITE(*,*)'  -Set max number of iterations for final'
       WRITE(*,*)'   fit (Default=',MITFIN,')'
       WRITE(*,*)'     Either set new value or 0 for default: '
           READ (THING,*) TOLS(4)
           IF (TOLS(4).EQ.0) THEN
             TOLS(4)=MITFIN
           ENDIF
       ELSEIF (TEMP.EQ.3) THEN
         WRITE (*,*)'What level of space efficiency'
         WRITE (*,*)'0. None'
         WRITE (*,*)'1. Moderate'
         WRITE (*,*)'2. Extreme: '
          READ (THING,*) FLAGS(11)
       ELSEIF (TEMP.EQ.4) THEN
         IF (FLAGS(8).EQ.1) THEN
          WRITE(*,*) 'Not applicable in current analysis'
         ELSE
         WRITE(*,*)'Do you want to output, to a separate file, the'
          WRITE (*,*) ' 0- nothing from this list'
          WRITE (*,*) ' 1- parameter estimates'
          WRITE (*,*) ' 2- point likelihoods: '
          READ (THING,*) FLAGS(20)
          IF (FLAGS(20).NE.0) THEN
           WRITE(*,*) 'What do you wish this file to be called: '
           READ (THING,'(A)') P3FYLE
           OPEN (UNIT=29,FILE=P3FYLE,STATUS = 'UNKNOWN')
          ENDIF
          WRITE (*,*) 'Do you want to output the indices and'
          WRITE (*,*) ' resulting allocations (0-no, 1=yes): '
          READ (THING,*) FLAGS(15)
          IF (FLAGS(15).NE.0) THEN
           WRITE(*,*) 'What do you wish this file to be called: '
           READ (THING,'(A)') P1FYLE
           OPEN (UNIT=56,FILE=P1FYLE,STATUS = 'UNKNOWN')
          ENDIF
          IF (FLAGS(10).NE.0) THEN
           WRITE (*,*) 'Do you want to output the bootstrap'
           WRITE (*,*) ' distribution values (0-no, 1-yes): '
           READ (THING,*) FLAGS(16)
           IF (FLAGS(16).NE.0) THEN
            WRITE(*,*) 'What do you wish this file to be called: '
            READ (THING,'(A)') P2FYLE
            OPEN (UNIT=57,FILE=P2FYLE,STATUS = 'UNKNOWN')
           ENDIF
          ENDIF
         ENDIF
         ELSEIF (TEMP.EQ.5) THEN
         IF (FLAGS(8).EQ.1) THEN
          WRITE(*,*) 'Not applicable in current analysis'
         ELSE
          FLAGS(12)=1-FLAGS(12)
           WRITE (*,*) ' (append entity group list to end of data file)'
         ENDIF
         ELSEIF (TEMP.EQ.6) THEN
          IF (FLAGS(17).GT.0) THEN
            FLAGS(17)=0
          ELSE
            FLAGS(17)=1
            IF (RANDTYPE.EQ.0) THEN
             WRITE (*,*) 'Unable to calculate Standard Errors as'
             WRITE (*,*) 'random number generator does not work'
             FLAGS(17)=0
            ELSE
             NEDRAN=1
             WRITE (*,*)'Which method of estimation'
             WRITE (*,*)' 0 Parametric'
             WRITE (*,*)' 1 Sampling with replacement'
             WRITE (*,*)' 2 weighted likelihood'
             IF (NCOV.EQ.2) THEN
               WRITE (*,*)' 4 information based method'
             ENDIF
             WRITE (*,*) 'Choose: '
             READ (THING,*) FLAGS(18)
             IF ((FLAGS(10).EQ.0).AND.(FLAGS(18).NE.4)) THEN
              WRITE (*,*) ' How many replications do you wish to use: '
              READ (THING,*) FLAGS(17)
             ELSE
              FLAGS(17)=FLAGS(10)
             ENDIF
             WRITE (*,*)' Enter name of output file for',
     &       ' Standard Errors: '
             READ (THING,'(A)') OFYLE2
             OPEN(UNIT=26,FILE=OFYLE2,STATUS='UNKNOWN')
            ENDIF
          ENDIF
          ELSEIF (TEMP.EQ.7) THEN
C         Option in New version
           WRITE (*,*)
           WRITE (*,*) 'This option is unavailable in this'
           WRITE (*,*) 'current version'
           WRITE (*,*)
          ELSEIF (TEMP.EQ.8) THEN
      IF (((FLAGS(4).EQ.2).OR.(NG1.EQ.1)).AND.(FLAGS(8).NE.1)) THEN
          FLAGS(6)=1-FLAGS(6)
         ELSE
          WRITE(*,*) 'Not applicable in current analysis'
         ENDIF
          ELSEIF (TEMP.EQ.9) THEN
         IF ((FLAGS(8).EQ.1).OR.(FLAGS(10).GT.0)) THEN
          WRITE(*,*) 'Not available when using bootstrap test'
         ELSE
           IF (FLAGS(7).EQ.0) THEN
      WRITE(*,*)'  1-Fixed user-defined degrees of freedom NU'
      WRITE(*,*)'    for each component'
      WRITE(*,*)'  2-Degrees of freedom NU estimated for each'
      WRITE(*,*)'    component (from user-supplied initial value)'
      WRITE(*,*)'  3-Common degrees of freedom NU estimated for the'
      WRITE(*,*)'    components (from user-supplied initial value)'
      WRITE(*,*)'  4-Degrees of freedom NU estimated for each component'
      WRITE(*,*)'    (moments estimates used as the initial values): '
             READ (THING,*) FLAGS(7)
             IF ((FLAGS(7).EQ.1).OR.(FLAGS(7).EQ.2)) THEN
             DO 533 II=1,NG1
532           WRITE (*,*) '   What value for NU',II,': '
              READ (THING,*) XUU(II)
              IF (XUU(II).LE.0) THEN
             WRITE(*,*)' ERROR: NU value must be greater than 0'
               GOTO 532
             ENDIF
533          CONTINUE
             ELSEIF (FLAGS(7).EQ.3) THEN
              WRITE (*,*) '   What value for NU: '
              READ (THING,*) XUU(1)
             DO 534 II=2,NG1
                XUU(II)=XUU(1)
534          CONTINUE
             ENDIF
           ELSE
             FLAGS(7)=0
           ENDIF
          ENDIF
          ELSEIF (TEMP.EQ.10) THEN
            IF (FLAGS(10).GT.0) THEN
            FLAGS(21)=1-FLAGS(21)
            ELSE
             WRITE (*,*) 'Only available when bootstrapping'
            ENDIF
          ELSEIF (TEMP.EQ.0) THEN
            GOTO 600
          ELSE
           WRITE(*,*)'Invalid menu number please select again'
          ENDIF
          GOTO 700

600       RETURN
          END

      SUBROUTINE USRPARAMETERS(NIND,NATT,NG,NCOV,XMU,XVAR,
     &           T,RNATT,RV)
C     This section is appropriate for FLAGS(4) = 2
C     where the user has supplied starting parameters for the EM
C     algorithm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,RNATT,RV(MNATT)
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION XMU(MAXNG,MNATT),XVAR(MAXNG,MNATT,MNATT),
     &          T(MAXNG)
      DO 400 K=1,NG
        READ (21,*) (XMU(K,J),J=1,NATT)
      WRITE (*,405) K
      WRITE (FYLENO,405) K
405   FORMAT (/2X,'User mean (as a row vector) for component ',I2)
      WRITE (*,415) (XMU(K,J),J=1,NATT)
      WRITE (FYLENO,415) (XMU(K,J),J=1,NATT)
415   FORMAT (2X,5G14.6)
        DO 404 I=1,NATT
         READ (21,*) (XVAR(K,I,J),J=1,I)
404     CONTINUE
400   CONTINUE
      READ (21,*) (T(K),K=1,NG)
C     Test if a common covariance matrix is specified (NCOV = 1)
      IF (NCOV.EQ.1) THEN
        WRITE (*,417)
        WRITE (FYLENO,417)
417     FORMAT (/2X,'User common component-covariance matrix ')
        DO 430 J=1,NATT
          DO 420 I=J+1,NATT
420         XVAR(1,J,I)=XVAR(1,I,J)
          WRITE (*,415) (XVAR(1,I,J),I=1,J)
          WRITE (FYLENO,415) (XVAR(1,I,J),I=1,J)
430     CONTINUE
      ELSE
        DO 450 K=1,NG
          WRITE (*,445) K
          WRITE (FYLENO,445) K
445       FORMAT (/2X,'User covariance matrix for component ',I2)
          DO 450 J=1,NATT
            DO 460 I=J+1,NATT
460           XVAR(K,J,I)=XVAR(K,I,J)
            WRITE (*,415) (XVAR(K,I,J),I=1,J)
            WRITE (FYLENO,415) (XVAR(K,I,J),I=1,J)
450     CONTINUE
      ENDIF
      WRITE (*,455)
      WRITE (FYLENO,455)
455   FORMAT (/2X,'User mixing proportion for each component')
      WRITE (*,465) (T(K),K=1,NG)
      WRITE (FYLENO,465) (T(K),K=1,NG)
465   FORMAT (5X,10F7.3)

C    Variable subset
      IF (FLAGS(19).EQ.2) THEN
        CALL REDUCE(NATT,NG,NCOV,XMU,XVAR,DV,V,RNATT,RV)
        WRITE(FYLENO,*)
        WRITE(FYLENO,*) ' Subset parameters'
        CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
      ENDIF

      RETURN
      END


      SUBROUTINE USRALLOC(NIND,NATT,NG,IDT,IER)
C     This subroutine sets up the initialisation for the EM algorithm
C     when an initial partition is given by the user (FLAGS(4) = 1)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      INTEGER IDT(MNIND)
      READ (21,*) (IDT(I),I=1,NIND)
      DO 300 I=1,NIND
       IF (IDT(I).GT.NG) THEN
         IER=6
c      ERROR as the input file may be in the wrong format of the
c      specified parameters do not match the input file.
         WRITE (*,*) 'READING ERROR!'
         WRITE (*,*) ' Are you sure the input file is in'
         WRITE (*,*) ' right form as ',IDT(I)
         WRITE (*,*) ' was read as a group number and you'
         WRITE (*,*) ' only have',NG,' groups.'
         WRITE (*,*) ' HINT: You might have typed the wrong'
         WRITE (*,*) '       value for number of points is'
         WRITE (*,*) '       ',NIND,' correct'

         RETURN
        ENDIF
300   CONTINUE
      IF ((FLAGS(18).NE.1).AND.(FLAGS(18).NE.2).AND.(FLAGS(3).NE.2))
     & THEN
      WRITE (FYLENO,305)
305   FORMAT (/2X,'Initial grouping as specified by input')
      WRITE (FYLENO,315) (IDT(I),I=1,NIND)
      WRITE (*,*) ' Partition read in = ',idt(1),idt(2),'...',idt(NIND)
315   FORMAT (2X,10I4)
      ENDIF

      RETURN
      END

C
C
C
      SUBROUTINE PRED(NIND,NATT,NCOV,NG,X,XMU,XVAR,V,DV,T,USA,WL,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      INTEGER USA(MNIND)
      COMMON /STORE2/ FLAGS,FYLENO
      EXTERNAL AUTOPARTITION
      EXTERNAL RANDNUM,DETERRANDOM
      DOUBLE PRECISION RANDNUM
      DIMENSION X(MNIND,MNATT),XVAR(MAXNG,MNATT,MNATT),T(MAXNG),
     &         DV(MAXNG),V(MAXNG,MNATT,MNATT),XMU(MAXNG,MNATT),
     &         IDT(MNIND),XUU(MAXNG),XLOGL(MITFIN),U(MNIND,MAXNG),
     &         W(MNIND,MAXNG),WL(MNIND),DEN(MNIND,MAXNG+1),
     &         XMAH(MNIND,MAXNG)
       IER=0
       IOUNT=1 
      CALL ESTEP(NIND,NATT,NG,X,XMU,V,T,DEN,WL,W,XUU,USA,DV,
     &            XLOGL,IOUNT,XMAH,IER)
      WRITE (FYLENO,*) ' Log Likelihood is ',XLOGL(IOUNT)
      IF (IER.EQ.-111) THEN
        WRITE(*,*) 'WARNING : Some points have zero Likelihood'
        WRITE(*,*) '         (will denote with 0 in grouping)'
        IER=0
      ENDIF
      CALL OUTLOOP(NIND,NATT,NG,XMU,DV,T,NCOV,IOUNT,XLOGL,
     &                   W,IDT,X,DEN,USA,U)

      RETURN
      END


      SUBROUTINE DISC(NIND,NATT,NCOV,NG,X,XMU,XVAR,V,DV,T,USA,WL,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      INTEGER USA(MNIND)
          COMMON /STORE2/ FLAGS,FYLENO
      EXTERNAL AUTOPARTITION
      EXTERNAL RANDNUM,DETERRANDOM
      DOUBLE PRECISION RANDNUM
      DIMENSION X(MNIND,MNATT),XVAR(MAXNG,MNATT,MNATT),T(MAXNG),
     &         DV(MAXNG),V(MAXNG,MNATT,MNATT),XMU(MAXNG,MNATT),
     &         IDT(MNIND),XUU(MAXNG),DEN(MNIND,MAXNG+1),
     &         W(MNIND,MAXNG),WL(MNIND),XLOGL(MITFIN),U(MNIND,MAXNG),
     &         XMAH(MNIND,MAXNG)
       IER=0
C     Initialisation of partition to zero for discrimination option
       DO 2 II=1,NIND
        IDT(II)=0
2      CONTINUE
      
C      Calculate parameter estimates
       CALL ESTIMATES(NIND,NATT,NG,X,IDT,WL,NCOV,XMU,V,XVAR,DV,
     &                  T,USA,IER)
       WRITE(FYLENO,*) '    Parameters from classified data'
       WRITE(FYLENO,*) '    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
C      Write parameter estimates to output file
       CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
      
       IF (IER.GT.0) THEN
         WRITE (FYLENO,*)'  Unable to continue covariance singular'
         WRITE (FYLENO,*)'  too few points in one of the components'
         RETURN
       ENDIF
      IOUNT=1
      CALL ESTEP(NIND,NATT,NG,X,XMU,V,T,DEN,WL,W,XUU,USA,DV,
     &            XLOGL,IOUNT,XMAH,IER)

      WRITE (FYLENO,*) ' Log Likelihood is ',XLOGL(IOUNT)
      IF (IER.EQ.-111) THEN
       WRITE(FYLENO,*)'  Warning: Some points have zero Likelihood'
       WRITE(FYLENO,*) '         (will denote with 0 in grouping)'
       WRITE(*,*) 'Warning : Some points have zero Likelihood'
      ENDIF
       CALL CAPART(NIND,NATT,NG,W,IDT,USA,XCC)
      WRITE(FYLENO,*)
      WRITE(FYLENO,*) '     *******************************'
      WRITE(FYLENO,*) '    FIT USING PARTIAL CLASSIFIED DATA'
      WRITE(FYLENO,*) '     *******************************'
      WRITE(FYLENO,*)'  Implied grouping for all unclassified entities'
      WRITE(FYLENO,*)' (with component membership of classified'
      WRITE(FYLENO,*)'  entities as specified)'
        WRITE(FYLENO,1177) (IDT(III),III=1,NIND)
        WRITE (FYLENO,*)
1177    FORMAT (2X,10I4)
        CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
        FLAGS(12)=0
        CALL ESTEP(NIND,NATT,NG,X,XMU,V,T,DEN,WL,W,XUU,USA,DV,
     &  XLOGL,IOUNT,XMAH,IER)
        IF (IER.GT.0) RETURN 
        XTMP=1
        CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,DV,
     &           XMU,WTOT,T,W,XUU,XMAH,XTMP,U,IER)
        IF (IER.GT.0) RETURN 
        FLAGS(12)=1
        WRITE(FYLENO,*)
        CALL OUTLOOP(NIND,NATT,NG,XMU,DV,T,NCOV,IOUNT,XLOGL,
     &               W,IDT,X,DEN,USA,U)

      RETURN
      END
C
C
C
      SUBROUTINE MAIN(NIND,NATT,NG0,NG1,NCOV,X,TOLS,USA,
     &   RV,RNATT,SIG,TXUU,TXMU,TXVAR,TT,TIDT,W,
     &   XMUS,XVARS,TS,XUU,AIC,BIC,AWE,TLL,PVAL,IDTS,SET,SEU,SEV,IER)

C   PURPOSE 
C     This is the main subroutine which controls the program 
C     and branches according to the options specified by FLAGS.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C    Constants that define array sizes at compilation time.
      INCLUDE 'EMMIX.max'

C    Global Parameters
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO

C   INPUT PARAMETERS
      INTEGER   FLAGS(40),     ! Flags summarising options chosen (see below)
     &          FYLENO,        ! File id for output file
     &          NIND,          ! Number of points (samples, or observations)
     &          NATT,          ! Number of dimensions (variates or attributes)
     &          NG0,           ! Minimum number of groups to be tested
     &          NG1,           ! Maximum number of groups to be tested
     &          NCOV,          ! Structure of covariance matrices
     &          IX,IY,IZ       ! Random seeds

      DIMENSION X(MNIND,MNATT) ! Data or sample. 

C   OTHER INPUT PARAMETERS
      INTEGER   TIDT(MNIND),              ! User defined partition of sample
     &          USA(MNIND),               ! Grouping of classified sample
     &          RNATT,                    ! Reduced number of dimensions
     &          RV(MNATT)                 ! Array of dimensions to be used
      DIMENSION TXUU(MAXNG),              ! User defined NU for fitting t-dist
     &          TXVAR(MAXNG,MNATT,MNATT), ! User defined covariance matrices
     &          TT(MAXNG),                ! User defined mixing proportions
     &          TXMU(MAXNG,MNATT),        ! User defined group means
     &          W(MNIND,MAXNG),           ! User defined posterior probabilities
     &           TOLS(4)                  ! User stopping tolerances for EM

C   OUTPUT PARAMETERS
      INTEGER    IDTS(MNIND,MAXNG)                !Partitions of the data
      DIMENSION  XMUS(MAXNG,MAXNG,MNATT),         !Group means
     &           XVARS(MAXNG,MAXNG,MNATT,MNATT),  !Group covariance matrices
     &           XUU(MAXNG),                      !Group NU values (for t-dist)
     &           TS(MAXNG,MAXNG),                 !Group mixing proportions
     &           AIC(MAXNG),                      !Akaike Information Criterion
     &           BIC(MAXNG),                      !Bayesian  "  "       "  "
     &           AWE(MAXNG),                      !Approx. Weight. Evidence
     &           PVAL(MAXNG),                     !P-value of -2log(Lambda)
     &           TLL(MAXNG)                       !-2log(Lambda)
C   OTHER OUTPUT PARAMETERS
      DIMENSION  SET(MAXNG),            !Standard Errors of mixing proportions
     &           SEU(MAXNG,MNATT),      !Standard Errors of group means
     &           SEV(MAXNG,MNATT,MNATT) !Standard Errors of covariance matrices
 

C   WORK PARAMETERS 
      INTEGER TFLAG
      DIMENSION XVAR(MAXNG,MNATT,MNATT),T(MAXNG),DV(MAXNG),
     &          XMU(MAXNG,MNATT),IDT(MNIND),V(MAXNG,MNATT,MNATT),
     &          TXML(MAXNG),XLOGL(MITFIN),FSEED(3),DEN(MNIND,MAXNG+1),
     &          XMAH(MNIND,MAXNG),WL(MNIND)
 
C     Store random seeds for the record 
      FSEED(1)=IX
      FSEED(2)=IY
      FSEED(3)=IZ

      IF (FLAGS(14).NE.1) THEN 
        DO 10 I=1,NIND
          WL(I)=1
10      CONTINUE
      ENDIF

C    Variable subset
      IF (FLAGS(19).EQ.2) THEN
        CALL REDUCE(NATT,NG0,NCOV,TXMU,TXVAR,DV,V,RNATT,RV)
        WRITE(FYLENO,*)
        WRITE(FYLENO,*) ' Subset parameters'
        CALL OUTESTIMATES(NIND,NATT,NG,NCOV,TXMU,V,TXVAR,DV,TT,TXUU)
      ENDIF

C    Invert user supplied covariance matrices
      IF (FLAGS(23).EQ.1) THEN
        CALL INVRT(NATT,NCOV,NG0,TXVAR,V,DV,IER,NULL)
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Generate a sample (OPTION 0)
      IF (FLAGS(8).EQ.0) THEN
        FYLENO=22
        CALL MVNG(NIND,NATT,NG0,TXMU,TXVAR,TT,X,IDT,IER)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Stand alone bootstrap or SE analysis (OPTION 1)
      ELSEIF (FLAGS(8).EQ.1) THEN
        FLAGS(3)=2
        FYLENO=22
C       If using a bootstrap based method 
        IF ((FLAGS(10).GT.0).OR.(FLAGS(18).NE.4)) THEN
          CALL MMOUT(NIND,NATT,NG0,X,TXMU,TXVAR,TT,NCOV,TXUU,TOLS,
     &    PVAL(1),TLL(1),USA,TIDT,SET,SEU,SEV,FSEED,RV,RNATT)
          FYLENO=22
          FLAGS(3)=1
        ENDIF
C       The Information-Based SE estimation
        IF (FLAGS(18).EQ.4) THEN
         ITEMP=FYLENO
           IER=0
          FYLENO=26
           CALL GDET(NATT,NG0,TXVAR,V,DV,IER,NULL)
           IF (IER.NE.0) THEN
            WRITE (FYLENO,*) '  Problem:'
            WRITE (FYLENO,*) '    Singular covariance matrix for'
            WRITE (FYLENO,*) '    user provided solution, unable'
            WRITE (FYLENO,*) '    to calculate standard errors.'
           ELSE
C        Calculate posterior probabilities W
         IOUNT=1
         CALL ESTEP(NIND,NATT,NG0,X,TXMU,V,TT,DEN,WL,W,TXUU,USA,DV,
     &            XLOGL,IOUNT,XMAH,IER)
C        Calculate standard errors via Information-based method
         CALL CALINFO(NIND,NATT,NG0,X,TXMU,V,W,TT,NCOV,TIDT,SET,SEU,SEV)
         ITEMP=FYLENO
C        Display Standard Errors
         CALL SEDISP(NATT,NG0,SEU,SEV,SET)
         FYELNO=ITEMP
         ENDIF
        ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Clustering for a single specified NG (OPTION 2)
ccJK        ELSEIF ((FLAGS(8).EQ.2).OR.(FLAGS(8).EQ.5)) THEN
        ELSEIF (FLAGS(8).EQ.2) THEN
        NG=NG0
        ING=1
        FYLENO=22
        FLAGS(3)=1
        DO 1815 II=1,NIND
          IDT(II)=TIDT(II)
1815    CONTINUE
        DO 1816 KK=1,NG
          T(KK)=TT(KK)
          XUU(KK)=TXUU(KK)
          DO 1816 II=1,NATT
            XMU(KK,II)=TXMU(KK,II)
            DO 1816 JJ=1,NATT
              XVAR(KK,II,JJ)=TXVAR(KK,II,JJ)
1816    CONTINUE

c      Call main clustering subroutine
          IER=0
          CALL NMM(NIND,NATT,NG,NCOV,IDT,W,X,WL,
     &     TXML(ING),DV,V,TOLS,XMU,XVAR,T,XUU,USA,FSEED,IER)
          IF (IER.EQ.6) RETURN 
c         Store results
          DO 3816 KK=1,NG
            TS(ING,KK)=T(KK)
            DO 3816 II=1,NATT
              XMUS(ING,KK,II)=XMU(KK,II)
              DO 3816 JJ=1,NATT
                XVARS(ING,KK,II,JJ)=XVAR(KK,II,JJ)
3816      CONTINUE
c      Estimate criteria values AIC,BIC etc.
          CALL CRITERIA(NG,TXML(ING),NIND,NATT,NCOV,AIC(ING),
     &                  BIC(ING),AWE(ING))
   	  WRITE (FYLENO,*) 'Criteria for this Clustering are'
	  WRITE (FYLENO,*) ' AIC  BIC'
	  WRITE (FYLENO,*) AIC(ING),BIC(ING)
          DO 1817 II=1,NIND
            IDTS(II,ING)=IDT(II)
1817       CONTINUE
         IF ((FLAGS(17).GT.0).AND.(FLAGS(18).NE.4)) THEN
c            call Bootstrapping SE subroutine
              CALL MMOUT(NIND,NATT,NG,X,XMU,XVAR,T,NCOV,XUU,TOLS,
     &   PVAL(ING),TLL(ING),USA,TIDT,SET,SEU,SEV,FSEED,RV,RNATT)
         ELSEIF (FLAGS(18).EQ.4) THEN
C          Calculate standard errors via Information-based method
           WRITE (*,*) 'Calculating standard errors...'
           ITEMP=FYLENO
           FYLENO=26
           IER=0
           CALL GDET(NATT,NG,XVAR,V,DV,IER,NULL)
           IF (IER.NE.0) THEN
            WRITE (FYLENO,*) '  Problem:'
            WRITE (FYLENO,*) '    Singular covariance matrix for'
            WRITE (FYLENO,*) '    solution provided by EMMIX, unable'
            WRITE (FYLENO,*) '    to calculate standard errors.'
           ELSE
           CALL CALINFO(NIND,NATT,NG0,X,XMU,V,W,T,NCOV,TIDT,SET,SEU,
     &                  SEV)
C        Display Standard Errors
         CALL SEDISP(NATT,NG0,SEU,SEV,SET)
           FYLENO=ITEMP
           ENDIF
         ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Discriminant rule with densities (OPTION 4)
      ELSEIF (FLAGS(8).EQ.4) THEN
       CALL DISC(NIND,NATT,NCOV,NG0,X,TXMU,TXVAR,V,DV,TT,USA,WL,IER)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Prediction (OPTION 5)
      ELSEIF (FLAGS(8).EQ.5) THEN
       CALL PRED(NIND,NATT,NCOV,NG0,X,TXMU,TXVAR,V,DV,TT,USA,WL,IER)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Produce estimates from partition (OPTION 6)
      ELSEIF (FLAGS(8).EQ.6) THEN
        NG=NG0
C       Calculate parameter estimates from allocation
        CALL ESTIMATES(NIND,NATT,NG,X,TIDT,WL,NCOV,XMU,V,XVAR,DV,
     &                  T,USA,IER)
        IOUNT=1
        CALL ESTEP(NIND,NATT,NG,X,XMU,V,T,DEN,WL,W,XUU,USA,DV,
     &  XLOGL,IOUNT,XMAH,IER)
        WRITE (FYLENO,*) 'Log_likelihood=',XLOGL(IOUNT)
C       Display parameter estimates to output file if required
         CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
         IF (IER.GT.0) THEN
           WRITE (FYLENO,*) 'ERROR: One or more of the resulting'
           WRITE (FYLENO,*) '       covariance matrices is singular'
           WRITE (*,*) 'ERROR: One or more of the resulting'
           WRITE (*,*) '       covariance matrices is singular'
           RETURN 
         ENDIF
         
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Fit a mixture model for a range of NG0 to NG1 (Option 3)
      ELSE  
          ING=0
        DO 9999 NG=NG0,NG1
          ING=ING+1
          FYLENO=22
          FLAGS(3)=1
          WRITE (FYLENO,*) '**************************************',
     &                     '***********************'
          WRITE (FYLENO,*) '    Results for ',NG,' cluster(s)'
          WRITE (FYLENO,*) '**************************************',
     &                     '***********************'
          WRITE (*,*) '-------------------------------------------',
     &                '--------'
          WRITE (*,11) NG
11        FORMAT (2X,'Starting Analysis for',I3,' cluster(s)')
          DO 815 II=1,NIND
            TIDT(II)=IDT(II)
815       CONTINUE
          DO 816 KK=1,NG
            XUU(KK)=TXUU(KK)
            TT(KK)=T(KK)
            DO 816 II=1,NATT
              TXMU(KK,II)=XMU(KK,II)
              DO 816 JJ=1,NATT
                TXVAR(KK,II,JJ)=XVAR(KK,II,JJ)
816       CONTINUE

c         Call main clustering subroutine
          IER=0
          CALL NMM(NIND,NATT,NG,NCOV,IDT,W,X,WL,
     &     TXML(ING),DV,V,TOLS,XMU,XVAR,T,XUU,USA,FSEED,IER)
          IF (IER.EQ.6) RETURN 
c         Store results
          DO 2816 KK=1,NG
            TS(ING,KK)=T(KK)
            DO 2816 II=1,NATT
              XMUS(ING,KK,II)=XMU(KK,II)
              DO 2816 JJ=1,NATT
                XVARS(ING,KK,II,JJ)=XVAR(KK,II,JJ)
2816       CONTINUE


c         Calculate various criteria AIC,BIC etc
          CALL CRITERIA(NG,TXML(ING),NIND,NATT,NCOV,AIC(ING),
     &                BIC(ING),AWE(ING))
   	  WRITE (FYLENO,*) 'Criteria for this Clustering are'
	  WRITE (FYLENO,*) ' AIC  BIC'
	  WRITE (FYLENO,*) AIC(ING),BIC(ING)

C         Record partition for the record
          DO 817 II=1,NIND
            IDTS(II,ING)=IDT(II)
817       CONTINUE
          
          IF (ING.GT.1) THEN
c          Calculate likelihood ratio statistic
            TLL(ING)=(-2)*(TXML(ING-1)-TXML(ING))

c          Approximate P-value of likelihood ratio by resampling
c          and/or estimate Standard errors of parameters
            IF ((FLAGS(10).GT.0).OR.(FLAGS(17).GT.0)) THEN
              FLAGS(3)=2
c            call Bootstrapping subroutine
              CALL MMOUT(NIND,NATT,NG-1,X,TXMU,TXVAR,TT,NCOV,XUU,TOLS,
     &   PVAL(ING),TLL(ING),USA,TIDT,SET,SEU,SEV,FSEED,RV,RNATT)
              FYLENO=22
	      WRITE (FYLENO,*) 'and P-value from bootstrap of',
     &                          PVAL(ING)
              FLAGS(3)=1
            ENDIF
          ENDIF
 
c       Display Standard Errors if needed
          IF (((FLAGS(17).GT.0).AND.(ING.GT.1)).OR.
     &         FLAGS(18).EQ.4) THEN
              CALL SEDISP(NATT,NG-1,SEU,SEV,SET)
          ENDIF

C         Estimate standard errors for final NG in the range since
C         bootstrap is not done for this fit.
          IF ((FLAGS(17).GT.0).AND.(NG.EQ.NG1)) THEN
            DO 915 II=1,NIND
             TIDT(II)=IDT(II)
915         CONTINUE
            DO 916 KK=1,NG
              TT(KK)=T(KK)
              DO 916 II=1,NATT
                TXMU(KK,II)=XMU(KK,II)
                DO 916 JJ=1,NATT
                  TXVAR(KK,II,JJ)=XVAR(KK,II,JJ)
916         CONTINUE
            FLAGS(3)=2
C           TFLAG=FLAGS(8)
            TFLAG=FLAGS(10)
            FLAGS(10)=0
C           FLAGS(8)=1
            CALL MMOUT(NIND,NATT,NG,X,TXMU,TXVAR,TT,NCOV,XUU,TOLS,
     &    PVAL(ING),TLL(ING),USA,TIDT,SET,SEU,SEV,FSEED,RV,RNATT)
            FLAGS(3)=1
C           FLAGS(8)=TFLAG
            FLAGS(10)=TFLAG
            FYLENO=22
c           Display standard errors
            CALL SEDISP(NATT,NG,SEU,SEV,SET)
           ENDIF

c         Set p-value for G-1 versus G to 0 for G=1
          IF (NG.EQ.1) THEN
            PVAL(1)=0
            TLL(1)=0
          ENDIF

c         Check if the P-value is still significant
          IF (SIG.GT.0) THEN
            IF (PVAL(ING).GT.SIG/100) THEN
              WRITE (FYLENO,*) '******************************'
              WRITE (FYLENO,*) 'Insignificant P-value reached'
              WRITE (FYLENO,*) ' ',ING-1,' vs ',ING 
              WRITE (FYLENO,*) '******************************'
              NG1=NG
              GOTO 99999
            ENDIF
          ENDIF

9999  CONTINUE
99999 CONTINUE
      ENDIF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      WRITE (*,*) '---------------------------------------------------'
      FYLENO=22

c     Display the a summary of the results to determine the number of groups 
      IF (FLAGS(8).EQ.3) THEN
        CALL ANASUM(NIND,NATT,NG0,NG1,NCOV,TXML,TLL,PVAL,AIC,BIC,AWE)
      ENDIF

c     If required write output files for external plots
      CALL EXPLOT(NIND,NATT,NG0,NG1,IDTS,X)

      RETURN
      END 

      SUBROUTINE ANASUM(NIND,NATT,NG0,NG1,NCOV,TXML,TLL,PVAL,AIC
     &                  ,BIC,AWE)
c     
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION TXML(MAXNG),TLL(MAXNG),PVAL(MAXNG),
     &           AIC(MAXNG),BIC(MAXNG),AWE(MAXNG)
      
      WRITE(FYLENO,*)
      WRITE(FYLENO,*) ' ANALYSIS SUMMARY'
      WRITE(FYLENO,*) ' ~~~~~~~~~~~~~~~~'
      WRITE(FYLENO,*)
      IF (FLAGS(10).EQ.0) THEN
        WRITE(FYLENO,*)'--------------------------------------------',
     &  '-----------------'
        WRITE(FYLENO,*)'| NG | Log Like | -2logLAM |    AIC   |'
     &  ,'    BIC   |'
        WRITE(FYLENO,*)'----------------------------------------------',
     &  '---------------'
      ELSE
        WRITE(FYLENO,*)'---------------------------------------------',
     &  '-----------------------------'
        WRITE(FYLENO,*)'| NG | Log Like | -2logLAM |    AIC   |'
     &  ,'    BIC   |  P-VAL  |'
        WRITE(FYLENO,*)'---------------------------------------------',
     &  '-----------------------------'
      ENDIF
      ING=0
      DO 59 NG=NG0,NG1
        ING=ING+1
C        CALL CRITERIA(NG,TXML(ING),NIND,NATT,NCOV,AIC(ING),BIC(ING),
C     &                AWE(ING))
        IF (FLAGS(10).EQ.0) THEN
          IF (TXML(ING).EQ.0) THEN
            WRITE(FYLENO,66) NG
66          FORMAT(' |',I3,' | NO VAL   | NO VAL   | NO VAL   |',
     &      ' NO VAL   | NO VAL   |')
          ELSE
            WRITE(FYLENO,58) NG,TXML(ING),TLL(ING),AIC(ING),BIC(ING)
          ENDIF
          WRITE(FYLENO,*)'----------------------------------------',
     &    '---------------------'
        ELSE
          IF (TXML(ING).EQ.0) THEN
            WRITE(FYLENO,67) NG
67          FORMAT(' |',I3,' | NO VAL   | NO VAL   | NO VAL   |',
     &      ' NO VAL   | NO VAL    |')
          ELSEIF (ING.EQ.1) THEN
            WRITE(FYLENO,68) NG,TXML(ING),AIC(ING),BIC(ING)
68          FORMAT(' |',I3,' |',F9.2,' |    -     |',F9.2,' |',
     &        F9.2,' |   -      |')
          ELSEIF (PVAL(ING).GE.0) THEN
            WRITE(FYLENO,58) NG,TXML(ING),TLL(ING),AIC(ING),BIC(ING),
     &      PVAL(ING)
          ELSE
            WRITE(FYLENO,56) NG,TXML(ING),TLL(ING),AIC(ING),BIC(ING),
     &      ((-1)*PVAL(ING))
          ENDIF
          WRITE(FYLENO,*)'----------------------------------------',
     &    '----------------------------------'
        ENDIF
57      FORMAT(' |',I3,' |',F9.2,' |',F9.2,' |',F9.2,' |',
     &        F9.2,' |',F9.2,' |')
58      FORMAT(' |',I3,' |',F9.2,' |',F9.2,' |',F9.2,' |',
     &        F9.2,' |  ',F9.2,' |')
56      FORMAT(' |',I3,' |',F9.2,' |',F9.2,' |',F9.2,' |',
     &        F9.2,' | <',F9.2,' |')
59      CONTINUE
      RETURN
      END


      SUBROUTINE EXPLOT(NIND,NATT,NG0,NG1,IDTS,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION X(MNIND,MNATT),IDTS(MNIND,MAXNG)

      IF (FLAGS(15).NE.0) THEN
        DO 96 II=1,NIND
          WRITE (56,98) II,(IDTS(II,ING),ING=1,NG1-NG0+1)
98        FORMAT(I5,20I3)
96      CONTINUE
      ENDIF
      IF (FLAGS(22).EQ.1) THEN
      DO 60 III=1,NIND
      WRITE (39,65) (X(III,JJJ),JJJ=1,NATT)
60     CONTINUE
65      FORMAT(10F10.5)
       CLOSE(39)
      ENDIF
      RETURN
      END


      SUBROUTINE SUMRY(NIND,NATT,NG0,NG1,NCOV,TOLS,NEDRAN,FI,
     &                 RV,RNATT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,RV(MNATT),FI,RNATT
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO
      COMMON /STORE3/ OFYLE,INFYLE,P1FYLE,P2FYLE,P3FYLE,OFYLE2
      DIMENSION TOLS(4)
      CHARACTER INFYLE*20, OFYLE*20,OFYLE2*20
      CHARACTER P1FYLE*20,P2FYLE*20,P3FYLE*20
      WRITE(FI,*)'      --------------------------------------------'
      IF (FLAGS(8).EQ.0) THEN 
         WRITE(FI,*)'        EMMIX - Generated sample from Mixture'
      ELSEIF (FLAGS(8).EQ.1) THEN
        IF ((FLAGS(10).GT.0).AND.(FLAGS(17).GT.0)) THEN
         WRITE(FI,*)'        EMMIX - Stand Alone Bootstrap',
     &   '-Standard Error Analysis'
        ELSEIF (FLAGS(10).GT.1) THEN
         WRITE(FI,*)'        EMMIX - Stand Alone Bootstrap Analysis'
        ELSE
         WRITE(FI,*)'        EMMIX - Stand Alone Standard Error',
     &    ' Analysis'
        ENDIF
      ELSEIF (FLAGS(8).EQ.2) THEN
      WRITE(FI,*)'       EMMIX - Mixture Analysis(Clustering) of' 
      WRITE(FI,*)'                 Sample in File ',INFYLE
      WRITE(FI,*)'                 for a specified range of g' 
      ELSEIF (FLAGS(8).EQ.3) THEN
      WRITE(FI,*)'       EMMIX - Mixture Analysis(Clustering) of'
      WRITE(FI,*)'                (Range of G) Sample in File ',INFYLE
      ELSEIF (FLAGS(8).EQ.4) THEN
      WRITE(FI,*)'       EMMIX - Discriminant Analysis of' 
      WRITE(FI,*)'                 Sample in File ',INFYLE
      ENDIF

      WRITE(FI,*)'      --------------------------------------------'

      WRITE(FI,*)'              input file:',INFYLE
      WRITE(FI,*)'              output file:',OFYLE
      WRITE(FI,*)'      --------------------------------------------'
      
      IF ((FLAGS(8).NE.4).AND.(FLAGS(8).NE.0)) THEN
      WRITE(FI,*)'      Initial fit Tolerance:',TOLS(1)
      WRITE(FI,580) TOLS(2)
      WRITE(FI,*)'      FINAL fit Tolerance:',TOLS(3)
      WRITE(FI,581) TOLS(4)
580   FORMAT (7X,'Initial fit Max Iterations: ',F6.0)
581   FORMAT (7X,'FINAL fit Max Iterations: ',F6.0)
      WRITE(FI,*)'      --------------------------------------------'
      ENDIF
      WRITE(FI,*)'             number of entities:',NIND
      WRITE(FI,*)'             number of variables:',NATT
      IF ((FLAGS(8).EQ.2).OR.(FLAGS(8).EQ.0)) THEN
       WRITE(FI,*)'             number of components:',NG0
      ELSE
       WRITE(FI,*)'        Range of number of components is:'
       WRITE(FI,*)'         ',NG0,' to ',NG1,' components'
      ENDIF
      IF ((FLAGS(19).EQ.1).OR.(FLAGS(19).EQ.-1)) THEN
        WRITE(FI,*)'        Subset of variables used in analysis'
        WRITE(FI,*)'       ',(RV(II),II=1,RNATT)
        WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(8).NE.0) THEN
      IF (NCOV.EQ.1) THEN
        WRITE(FI,*)'                Equal covariance matrices'
      ELSEIF (NCOV.EQ.2) THEN
        WRITE(FI,*)'            Unrestricted covariance matrices'
      ELSEIF (NCOV.EQ.3) THEN
        WRITE(FI,*)'            Equal diagonal covariance matrices'
      ELSEIF (NCOV.EQ.4) THEN
        WRITE(FI,*)'        Unrestricted diagonal covariance matrices'
      ENDIF
      WRITE(FI,*) '      --------------------------------------------'
      IF (FLAGS(4).EQ.3) THEN
      WRITE(FI,*) '      automatic methods used for initial groupings'
      WRITE(FI,*)'      with ',FLAGS(26),
     &' random start(s), ',INT(FLAGS(5))
      WRITE(FI,*) '      k-mean(s) starts and'
      ELSEIF (FLAGS(4).EQ.2) THEN
       WRITE(FI,*) '      user-defined parameters used to start the' 
       WRITE(FI,*) '      EM algorithm' 
      ELSEIF (FLAGS(4).EQ.1) THEN
       WRITE(FI,*) '      user-defined grouping of sample used to'
       WRITE(FI,*) '      start the EM algorithm' 
      ELSEIF (FLAGS(4).EQ.4) THEN
       WRITE(FI,*) '      user-defined posterior probabilities used'
       WRITE(FI,*) '      to start the EM algorithm' 
      ENDIF
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF

      IF (NEDRAN.GT.0) THEN
       IF (RANDTYPE.EQ.1) THEN
         WRITE(FI,536) IX,IY,IZ
536      FORMAT(7X,'random seeds = ',I5,I5,I5)
       ELSE
         WRITE(FI,*) '      random seed = ',SEED
       ENDIF
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(10).GT.0) THEN
        WRITE(FI,*) '      Stepwise bootstrap analysis to be done  '
        WRITE(FI,*) '      with ',FLAGS(10),' replications'
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(2).EQ.1) THEN
        WRITE(FI,*) '      Stochastic EM algorithm used'
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(11).EQ.1) THEN
        WRITE(FI,*) '      Partial space saving operating'
      WRITE(FI,*) '      --------------------------------------------'
      ELSEIF (FLAGS(11).EQ.2) THEN
        WRITE(FI,*) '      Extreme space saving operating'
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(12).EQ.1) THEN
        WRITE(FI,*) '      Partial user grouping utilized'
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(7).EQ.1) THEN
        WRITE(FI,*) '      Mixture of : Ts user defined NUs'
      WRITE(FI,*) '      --------------------------------------------'
      ELSEIF (FLAGS(7).EQ.2) THEN
        WRITE(FI,*)'     Mixture of : Ts Iterative NU started at'
        WRITE(FI,*)'     user estimate with unequal NUs for all groups' 
      WRITE(FI,*) '      --------------------------------------------'
      ELSEIF (FLAGS(7).EQ.3) THEN
        WRITE(FI,*)'     Mixture of : Ts Iterative NU started at'
        WRITE(FI,*)'     user estimate with equal NUs for all groups' 
      WRITE(FI,*) '      --------------------------------------------'
      ELSEIF (FLAGS(7).EQ.4) THEN
        WRITE(FI,*)'     Mixture of : Ts Iterative NU started at'
      WRITE(FI,*)'     moment estimate with unequal NUs for all groups' 
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(15).EQ.1) THEN
        WRITE(FI,*) '      Plotting information of cluster analysis'
        WRITE(FI,*) '      sent to file -',P1FYLE
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(16).EQ.1) THEN
        WRITE(FI,*) '      Plotting information of bootstrap analysis'
        WRITE(FI,*) '      sent to file -',P2FYLE
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      IF (FLAGS(18).GT.0) THEN
        WRITE(FI,*) '      Standard Errors to be estimated and reported' 
      WRITE(FI,*) '      --------------------------------------------'
      ENDIF
      WRITE(FI,*)

597   CONTINUE

599    RETURN
      END







      SUBROUTINE MMOUT(NIND,NATT,NG,X,XMU,XVAR,T,NCOV,XUU,TOLS,
     &         PVAL,XVAL,USA,TIDT,SET,SEU,SEV,FSEED,RV,RNATT)
c      This subroutine displays all the relevant information from the
c      EM algorithm to the be used by the program MMresamp.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      EXTERNAL RESAMP
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM
      REAL LLAM
      INTEGER FLAGS(40),FYLENO,RV(MNATT),FI,USA(MNIND),RNATT,
     &        TIDT(MNIND)
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO
      COMMON /STORE3/ OFYLE,INFYLE,P1FYLE,P2FYLE,P3FYLE,OFYLE2
      DIMENSION T(MAXNG),XMU(MAXNG,MNATT),XVAR(MAXNG,MNATT,MNATT),
     &  IRANK(MAXREP),LLAM(MAXREP),TXML(2,MAXREP),SEEDS(MAXREP,3)
      DIMENSION IERS(MAXREP,2),DVS(MAXREP,MAXNG),TOLS(4),FSEED(3),
     &  X(MNIND,MNATT),XUU(MAXNG),
     &  SET(MAXNG),SEU(MAXNG,MNATT),SEV(MAXNG,MNATT,MNATT)
      CHARACTER INFYLE*20, OFYLE*20, OFYLE2*20
      CHARACTER P1FYLE*20,P2FYLE*20,P3FYLE*20
      CHARACTER OUTFLE*20
      
      WRITE (*,*) 'Running Bootstrap'
      WRITE (*,*) ' --------------------------------',
     & '----------------'
      IF ((FLAGS(8).NE.1).AND.(FLAGS(10).GT.0)) THEN
        I100=INT(FLOAT(NG)/100.0)
        I10=INT(FLOAT(NG-(I100*100))/10.0)
        I1=NG-I100*100-I10*10
        J100=INT(FLOAT(NG+1)/100.0)
        J10=INT(FLOAT(NG+1-(J100*100))/10.0)
        J1=NG+1-J100*100-J10*10
        IF (J100.GT.0) THEN
       OUTFLE='boot'//CHAR(I100+48)//CHAR(I10+48)//CHAR(I1+48)//
     &    'vs'//CHAR(J100+48)//CHAR(J10+48)//CHAR(J1+48)//'.out'
        ELSEIF (J10.GT.0) THEN
          OUTFLE='boot'//CHAR(I10+48)//CHAR(I1+48)//
     &    'vs'//CHAR(J10+48)//CHAR(J1+48)//'.out'
        ELSE
          OUTFLE='boot'//CHAR(I1+48)//
     &    'vs'//CHAR(J1+48)//'.out'
        ENDIF
           
        OPEN(UNIT=25,FILE=OUTFLE,STATUS='UNKNOWN')
      ENDIF  
C      SEEDS(1,1)=RANDNUM()*30000
C      SEEDS(1,2)=RANDNUM()*30000
C      SEEDS(1,3)=RANDNUM()*30000
       SEEDS(1,1)=FSEED(1)    
       SEEDS(1,2)=FSEED(2)    
       SEEDS(1,3)=FSEED(3)    
       IX=FSEED(1)
       IY=FSEED(2)    
       IZ=FSEED(3)    
      IF (FLAGS(17).GT.0) THEN
        NREP=FLAGS(17)
      ELSE
        NREP=FLAGS(10)
      ENDIF
      CALL RESAMP(NIND,NATT,NG,NCOV,X,XMU,XVAR,DVS,
     & T,TIDT,NREP,SEEDS,LLAM,TXML,TOLS,USA,
     & SET,SEU,SEV,FSEED,XUU,IERS,IER)
      IF (IER.EQ.12) THEN
       WRITE(*,*) 'FATAL ERROR: The random number generator does'
        WRITE(*,*) '             not seem to be working.'
        WRITE(*,*) ' Please modify or replace random number generator'
        WRITE(25,*) 'FATAL ERROR: The random number generator does'
        WRITE(25,*) '             not seem to be working.'
        WRITE(25,*) ' Please modify or replace random number generator'
      ENDIF

        IF ((FLAGS(17).GT.0).AND.(FLAGS(18).NE.4)) THEN
        IF (NREP.GT.1) THEN
            ITEMP=FYLENO
            FYLENO=26
              CALL SEDISP(NATT,NG,SEU,SEV,SET)
            FYLENO=ITEMP
        ELSE
         WRITE (*,*) 'Cannot calculate SE with only 1 iteration'
         WRITE (26,*) 'Cannot calculate SE with only 1 iteration'
        ENDIF
      ENDIF
      IF (FLAGS(10).EQ.0) RETURN

      IF (NREP.GT.1) THEN
        IF (FLAGS(8).EQ.1) THEN
          REWIND(25) 
        ELSE
          CLOSE(25)
          OPEN(UNIT=25,FILE=OUTFLE,STATUS='UNKNOWN')
        ENDIF
        IF (FLAGS(10).GT.0) THEN
          FI=25
          CALL SUMRY(NIND,NATT,NG,NG+1,NCOV,TOLS,1,FI,RV,RNATT)
        ENDIF
C       Sort the results into descending order
C       Call Sorting Subroutine
        WRITE(*,*) '  sorting ...'
        CALL SORT(NREP,LLAM,IRANK)
C       Write output to output file
        WRITE (25,55)
55      FORMAT('INDEX',7X,'LIK_H0',8X,'LIK_H1',7X,
     &  '-2logLAMDA',6X,'   %')
        PVAL=1/(REAL(NREP)+1)
        PERC=(-1)*PVAL
          IF ((FLAGS(11).EQ.2).AND.(FLAGS(8).NE.1)) THEN
            WRITE (25,*) '****************************************'
            WRITE (25,*) '* No Bootstrap replication Information *'
            WRITE (25,*) '* is outputed since the extreme space  *'
            WRITE (25,*) '*  saving option was chosen by use.    *'
            WRITE (25,*) '****************************************'
         ENDIF
        DO 125 I=1,NREP
          PERC=I/(REAL(NREP)+1)
          J=IRANK(I)
C         Find P-value of given statistic value
          IF (XVAL.LE.LLAM(J)) THEN
            PVAL=PERC
          ENDIF
          IF ((FLAGS(11).NE.2).OR.(FLAGS(8).EQ.1)) THEN
            IF ((IERS(J,1).EQ.0).AND.(LLAM(J).GT.0)) THEN
             WRITE (25,90) J,TXML(1,J),TXML(2,J),LLAM(J),PERC*100
90           FORMAT(I4,F15.3,F15.3,F15.3,F15.3,'%')
            ELSEIF (IERS(J,2).EQ.23) THEN
             WRITE(25,95) J,IERS(J,2)
95           FORMAT(I4,4X,
     &       'ERR: Problem in generation of sample (code',I2,')')
            ELSEIF (IERS(J,1).EQ.5) THEN
             WRITE(25,100) J,IERS(J,2)
100          FORMAT(I4,4X,
     &       'ERR: Problem with fit under H0 (code',I2,')')
            ELSEIF (IERS(J,1).EQ.6) THEN
             WRITE(25,105) J,TXML(1,J),IERS(J,2)
105          FORMAT(I4,F15.3,1X,
     &       'ERR:Auto start unable to find a start(code',I2,')')
            ELSEIF (IERS(J,1).EQ.7) THEN
             WRITE(25,110) J,TXML(1,J),IERS(J,2)
110          FORMAT(I4,F15.3,4X,
     &       'ERR: Problem with fit under H1(code',I2,')')
            ELSEIF (LLAM(J).LE.0) THEN
             WRITE (25,115) J,TXML(1,J),TXML(2,J),LLAM(J)
115          FORMAT(I4,F15.3,F15.3,F15.3,
     &       ' < Problem (LH0>LH1)')
            ENDIF
            WRITE (25,120) SEEDS(J,1),SEEDS(J,2),SEEDS(J,3)
120         FORMAT(4X,'Seeds= ',F15.1,F15.1,F15.1)
            WRITE (25,*) 'DET= ',(DVS(J,KK),KK=1,NG+1)
           ENDIF
125     CONTINUE
        CLOSE(25)
        IF (FLAGS(16).GT.0) THEN
          DO 130 I=1,NREP
            J=IRANK(I)
            WRITE (57,135) J,TXML(1,J),TXML(2,J),LLAM(J)
135         FORMAT(I4,F15.3,F15.3,F15.3)
130       CONTINUE      
        ENDIF
      ENDIF
      WRITE(*,*)'Bootstrap Complete'
      RETURN
      END

C     The following auxiliary subroutines are used

      SUBROUTINE SORT(N,IN,IRANK)
C     This subroutine, written by W.Whiten sorts a given array
C     It is used to sort the replicated values of -2log(Lambda)

      IMPLICIT REAL (A-H,O-Z)
       INCLUDE 'EMMIX.max'

      REAL IN
      INTEGER S
      DIMENSION IN(MAXREP),IRANK(MAXREP)
C     Set up index
      DO 200 I=1,N
         IRANK(I)=I
200   CONTINUE
C     Find initial step size
      STEP=1
      DO 210 S=1,N
        STEP=STEP*2
        IF (STEP.GT.N) GOTO 220
210   CONTINUE
220     STEP=(STEP-1)/2

C     Sort subsections until done with step of one
      DO 260 S=1,N
        STEP=STEP/2
        IF ((STEP.LT.1).AND.(STEP.GT.(0.6))) THEN
          STEP=1
        ELSEIF (STEP.LT.(0.6)) THEN
          GOTO 299
        ENDIF
C       Do step separate subsection sorts
        DO 250 I=1,STEP
C         Bubble sort of items step apart
          DO 240 J=I,(N-STEP+.5),STEP
            JS=J+STEP+.5
            IF (IN(IRANK(J)).LT.(IN(IRANK(JS)))) THEN
C             Not in order so move up
              IT=IRANK(JS)
              K=J
              KS=JS
              DO 230 M=1,N
                IRANK(KS)=IRANK(K)
                KS=K
                K=K-STEP+.5
                IF ((K.LT.1).OR.(IN(IRANK(K)).GT.IN(IT))) THEN
                  IRANK(KS)=IT
                  GOTO 240
                ENDIF
230           CONTINUE
            ENDIF
240       CONTINUE
245       CONTINUE
250       CONTINUE
260     CONTINUE
299   RETURN
      END 


      SUBROUTINE SEDISP(NATT,NG,SEU,SEV,SET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION SEU(MAXNG,MNATT),SEV(MAXNG,MNATT,MNATT),SET(MAXNG)
c       Display Standard Errors
          IF ((FLAGS(17).EQ.1).AND.(FLAGS(18).NE.4)) THEN
           WRITE (FYLENO,*) '  Unable to estimates SE'
           WRITE (FYLENO,*) '  from one replication'
          ELSE
           WRITE (FYLENO,*)
           WRITE (FYLENO,*) '           STANDARD ERRORS'
           WRITE (FYLENO,*)'           for ',NG,' groups'
         DO 917 K=1,NG
           WRITE (FYLENO,*)
           WRITE (FYLENO,*) '  SE of Mean for Cluster ',K
           WRITE (FYLENO,918) (SEU(K,J),J=1,NATT)
           WRITE (FYLENO,*) '  SE of Covariance for Cluster ',K
         DO 917 J=1,NATT
           WRITE (FYLENO,918) (SEV(K,J,L),L=1,J)
917      CONTINUE
           IF (NG.GE.2) THEN
            WRITE (FYLENO,*)
            WRITE (FYLENO,*) '  SE of Mixing Proportions'
            WRITE (FYLENO,918) (SET(K),K=1,NG-1)
           ENDIF
          ENDIF
918        FORMAT(2X,8F10.6)
         
          RETURN
          END
   
        SUBROUTINE INPSUM(NIND,NATT,NG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO
      COMMON /STORE3/ OFYLE,INFYLE,P1FYLE,P2FYLE,P3FYLE,OFYLE2
      CHARACTER INFYLE*20, OFYLE*20,OFYLE2*20
      CHARACTER P1FYLE*20,P2FYLE*20,P3FYLE*20
      
      WRITE(*,*)'       The Input File should be of the form:'
      IF ((FLAGS(8).GT.1).OR.(FLAGS(18).GE.1)) THEN
       WRITE(*,*)'       The sample consisting of'
       WRITE(*,*)'       ',NIND,' rows and',NATT,' columns'
      ENDIF
      IF (FLAGS(12).GT.0) THEN
       WRITE(*,*)'       Point group'
       WRITE(*,*)'       repeated once per line as required'
       WRITE(*,*)'       finish with -1 -1'
      ENDIF
      IF ((FLAGS(8).LE.1).OR.(FLAGS(4).EQ.2)) THEN
        WRITE(*,*)'       Mean component 1'
        WRITE(*,*)'       Covariance component 1'
        IF (NG.GT.1) THEN 
        WRITE(*,*)'       repeat for other',NG-1,' components'
        WRITE(*,*)'       ',NG,' mixing proportions'
        ELSE
        WRITE(*,*)'       a single 1'
        ENDIF
      ENDIF
      IF ((FLAGS(4).EQ.1).OR.(FLAGS(18).GE.1)) THEN
        WRITE(*,*)'       Allocations for each of',NIND,' points'
        WRITE(*,*)'         eg 1 1 2 2 1 2 1 2'
      ELSEIF (FLAGS(4).EQ.4) THEN
        WRITE(*,*)'       ',
     &           NG,'Posterior probabilities on each line' 
        WRITE(*,*) '      for',NIND,'points'
      ENDIF 
      WRITE(*,*)'      --------------------------------------------'
      RETURN
      END

      SUBROUTINE REDUCE(NATT,NG,NCOV,XMU,XVAR,DV,V,RNATT,RV)
C      This subroutine is used when a subset of the variables are
C      used for clustering but the parameter estimates read in
C      from the input file are for all the variables.
C      The unused variables are simply removed from the parameter
C      estimates.

       implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,RNATT,RV(MNATT)
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION XMU(MAXNG,MNATT),TXMU(MAXNG,MNATT),
     &    XVAR(MAXNG,MNATT,MNATT),TXVAR(MAXNG,MNATT,MNATT),
     &    DV(MAXNG)
        DO 1050 K=1,NG
          DO 1050 I=1,RNATT
           TXMU(K,I)=XMU(K,RV(I))
           IF ((NCOV.EQ.2).OR.(NCOV.EQ.0)) THEN
             DO 1048 J=1,RNATT
               TXVAR(K,I,J)=XVAR(K,RV(I),RV(J))
1048         CONTINUE
           ELSE
             DO 1049 J=1,RNATT
               TXVAR(K,I,J)=XVAR(1,RV(I),RV(J))
1049         CONTINUE
           ENDIF
1050   CONTINUE
        DO 1060 K=1,NG
          DO 1060 I=1,RNATT
           XMU(K,I)=TXMU(K,I)
             DO 1060 J=1,RNATT
               XVAR(K,I,J)=TXVAR(K,I,J)
1060   CONTINUE
       NATT=RNATT
       RETURN
       END

      SUBROUTINE INVRT(NATT,NCOV,NG,XVAR,V,DV,IER,NULL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION XVAR(MAXNG,MNATT,MNATT),DV(MAXNG),
     &  V(MAXNG,MNATT,MNATT)
 
         IF (NCOV.EQ.1) THEN
          CALL GDET(NATT,NG,XVAR,V,DV,IER,NULL)
          DO 4401 K=2,NG
           DV(K)=DV(1)
           DO 4401 J=1,NATT
            DO 4401 I=1,NATT
              XVAR(K,I,J)=XVAR(1,I,J)
              V(K,I,J)=V(1,I,J)
4401     CONTINUE

        ELSE
         CALL GDET(NATT,NG,XVAR,V,DV,IER,NULL)
        ENDIF
         IF (IER.NE.0) THEN
          WRITE (FYLENO,*) '  Problem:'
          WRITE (FYLENO,*) '    User defined covariance matrix'
          WRITE (FYLENO,*) '    is singular.'
         ENDIF
       RETURN
       END

C
C     
C     This group of subroutines handles the generation of initial partitions
C     for the EM algorithm.  This is done via various clustering methods plus
C     random starts
C     Written by D.Peel Oct 1994 based on program pkmmgen 1992

      SUBROUTINE AUTOPARTITION(NIND,NATT,NG,NCOV,X,
     &                  WBEST,WL,TOLS,XUU,USA,FSEED,IER)
C      This subroutine is the main subroutine and controls the calling and
C      interaction of the other subroutines to generate initial partitions

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
       INTEGER FLAGS(40),FYLENO,USA(MNIND)
       COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
       COMMON /STORE2/ FLAGS,FYLENO
       DIMENSION X(MNIND,MNATT),WBEST(MNIND,MAXNG),XML(MSTART),
     &           TOLS(4),WL(MNIND),FSEED(3),
     &           XUU(MAXNG),XUUBEST(MAXNG),CORIND(MSTART)
         CHARACTER*33 NAMBEST
        DOUBLE PRECISION RANDNUM
	EXTERNAL RANDNUM,KMEANS,HIER
	MAXNS=MSTART
        NS=0
c	WRITE (FYLENO,109)
	XMLBEST=0
	DO 12 I=1,MAXNS
	 CORIND(I)=0
12      CONTINUE	  
	FFG=-1
        WRITE(FYLENO,*) '     Results of search for initial Grouping'
        WRITE(FYLENO,*) '     ******* ** ****** *** ******* ********'
15       FORMAT (2X,'------------------------------------'
     &                ,'------------------------------------')
        IF (HIRFLG.EQ.1) THEN
         CALL HIERCON(NIND,NATT,NG,NCOV,X,WL,WBEST,XMLBEST,XML,NS,
     &        MAXNS,TOLS,NAMBEST,XUU,XUUBEST,USA,CORIND,FFG,IER)
         IF (FLAGS(11).EQ.0) THEN
         WRITE(FYLENO,*)
         WRITE(FYLENO,15)
         ENDIF
        ENDIF
        IF (FLAGS(2).EQ.1) THEN
         IX=FSEED(1)
         IY=FSEED(2)
         IZ=FSEED(3)
        ENDIF
        CALL RANDST(NIND,NATT,NG,NCOV,X,WL,WBEST,XMLBEST,XML,
     &    NS,MAXNS,TOLS,NAMBEST,XUU,XUUBEST,USA,CORIND,FFG,IER)
        IF (FLAGS(11).EQ.0) WRITE(FYLENO,15)
        CALL KMEANCON(NIND,NATT,NG,NCOV,X,WL,WBEST,XMLBEST,XML,NS,
     &        MAXNS,TOLS,NAMBEST,XUU,XUUBEST,USA,CORIND,FFG,IER)
        WRITE(FYLENO,15)
        
        IF (FFG.EQ.-1) THEN
	   IER=31
           IF (FLAGS(3).EQ.1) THEN
  	    WRITE(*,*)'  ERROR : Auto Start unable to find start'
           ENDIF
	   WRITE(FYLENO,*)'  ERROR :Auto Start unable to find start'
        ELSE
C          Output the likelihood values
           IER=0
           WRITE(FYLENO,25)
25        FORMAT(//,2X,'Final log likelihood values from each',
     &                  ' initial grouping'/)
	   IF (NS.GT.MAXNS-1) THEN
  	     WRITE (FYLENO,35) MAXNS
35	 FORMAT('Due to compilation restraints only the values from',I3,
     &	' starting grouping',/,'can be stored and listed. All',
     &  'values have been considered and this',/,' restriction',
     &  ' does not effect the programs performance in any way.',/,
     & 'To list all values increase matrix size and re-compile program')
           ENDIF
           WRITE(FYLENO,45)(XML(IS),IS=1,NS)
45         FORMAT(2X,6G12.4)
      WRITE(FYLENO,*)'  Best initial grouping (corresponding to the'
      WRITE(FYLENO,*)'  highest value of likelihood) found by the '
      WRITE(FYLENO,*)'  ',NAMBEST,' method'
           WRITE(FYLENO,15)
       ENDIF
       DO 50 KK=1,NG
        XUU(KK)=XUUBEST(KK)
50     CONTINUE
       RETURN
       END


       SUBROUTINE RANDST(NIND,NATT,NG,NCOV,X,WL,WBEST,XMLBEST,
     &     XML,NS,MAXNS,TOLS,NAMBEST,XUU,XUUBEST,USA,CORIND,FFG,IER)
C      This subroutine controls the generation of random starts and for
C      each random partition calls the subroutine FIT.

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       EXTERNAL RANDNUM
       DOUBLE PRECISION RANDNUM,R
       INCLUDE 'EMMIX.max'
       INTEGER FLAGS(40),FYLENO,USA(MNIND)
       COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
       COMMON /STORE2/ FLAGS,FYLENO
         DIMENSION X(MNIND,MNATT),XML(MSTART),TOLS(4),
     &     IDT(MNIND),WBEST(MNIND,MAXNG),W(MNIND,MAXNG),WL(MNIND),
     &     XUU(MAXNG),XUUBEST(MAXNG),CORIND(MSTART)
         CHARACTER*33 NAMBEST
         IF (FLAGS(11).EQ.0) THEN
         WRITE(FYLENO,105)
105      FORMAT(//2X,'           Random Initial Grouping ',/
     &           ,2X,'           ~~~~~~ ~~~~~~~ ~~~~~~~~')
         WRITE(FYLENO,*)  FLAGS(26),' Random Starts used '
         IF (FLAGS(1).EQ.100) THEN
         WRITE(FYLENO,*)  'All the data used to form starts'
         ELSE
      WRITE(FYLENO,*) FLAGS(1),' percent of the data used'
         ENDIF
         ENDIF
         IF (FLAGS(3).EQ.1) THEN
           WRITE(*,*)'Examining random starts...'
           WRITE(*,*)'Random start Number',
     &     '  ( Random seeds )'
           WRITE(FYLENO,*) 'Random start Number',
     &     '  ( Random seeds )         Like'
         ENDIF
         DO 199 JJ=1,FLAGS(26)
           IF (FLAGS(3).EQ.1) THEN
             IF (RANDTYPE.EQ.1.0) THEN
    	       WRITE(*,*) '   random start ',JJ,IX,IY,IZ
             ELSE
  	       WRITE(*,*) '   random start ',JJ,SEED
             ENDIF
           ENDIF
           IF (FLAGS(11).EQ.0) WRITE (FYLENO,175)
175       FORMAT (2X,'------------------------------------'
     &              ,'------------------------------------')
             IF (RANDTYPE.EQ.1.0) THEN
    	       WRITE(FYLENO,*) '  random start ',JJ,IX,IY,IZ
             ELSE
  	       WRITE(FYLENO,*) '  random start ',JJ,SEED
             ENDIF
cc          WRITE(FYLENO,*) '          Random start number',JJ
cc          WRITE(FYLENO,185)(IDT(MM),MM=1,NIND)
cc185       FORMAT(2X,10I4)
cc          IF (RANDTYPE.EQ.1) THEN
cc        WRITE(FYLENO,*) '  Random seeds for this start = ',IX,IY,IZ
cc          ELSE
cc            WRITE(FYLENO,*)'  Random seed for this start = ',IX
cc          ENDIF
cc          WRITE(FYLENO,*)
	   IF (FLAGS(1).EQ.100) THEN
             DO 110 I=1,NIND
	       R=RANDNUM()
               G=R*FLOAT(NG)
               JG=INT(G)+1
               IF (JG.GT.NG) JG=NG
               IDT(I)=JG
110         CONTINUE
          ELSE
             DO 111 I=1,NIND
               R=RANDNUM()
               XCUT=FLOAT(FLAGS(1))/100.00
               IF (R.GE.XCUT) THEN
                IDT(I)=0 
               ELSE
	        R=RANDNUM()
                G=R*FLOAT(NG)
                JG=INT(G)+1
                IF (JG.GT.NG) JG=NG
                IDT(I)=JG
               ENDIF
111         CONTINUE
          ENDIF
C           Now we have the data and an IDT initial allocation
190          NS=NS+1
C           Fit a multivariate normal mixture model to the data set X via
C           the EM algorithm initialised with the partition from random
C           assignment IDT
c          FLAGS(9)=3
          CALL FIT(NIND,NATT,NG,NCOV,X,WL,XUU,XUUBEST,USA,
     &	     IDT,WBEST,XML,XMLBEST,NS,MAXNS,TOLS,W,CORIND,FFG,IER)
   	  WRITE(FYLENO,685) XML(NS)
685     FORMAT(2X,'Log Likelihood value from EM algorithm started',
     &        /2X,'from this grouping is',F15.3)
        IF (XML(NS).EQ.XMLBEST) THEN
           I100=INT(FLOAT(JJ)/100.0)
           I10=INT(FLOAT(JJ-(I100*100))/10.0)
           I1=JJ-I100*100-I10*10
           NAMBEST='Random Start  '//CHAR(I100+48)
     &                            //CHAR(I10+48)
     &                            //CHAR(I1+48)
        ENDIF
199     CONTINUE
        RETURN
       END


       SUBROUTINE HIERCON(NIND,NATT,NG,NCOV,X,WL,WBEST,XMLBEST,
     &    XML,NS,MAXNS,TOLS,NAMBEST,XUU,XUUBEST,USA,CORIND,FFG,IER)
C     This subroutine calls the various Hierarchical clustering methods used
C     and for each resulting grouping calls the subroutine FIT.
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       EXTERNAL HIER
       INCLUDE 'EMMIX.max'
       INTEGER FLAGS(40),FYLENO,USA(MNIND)
       COMMON /STORE2/ FLAGS,FYLENO
         DIMENSION X(MNIND,MNATT),XML(MSTART),IDT(MNIND),
     &     WBEST(MNIND,MAXNG),BETA(MHIER),ISU(MHIER),IS(MHIER),
     &     TOLS(4),W(MNIND,MAXNG),WL(MNIND),XUU(MAXNG),
     &     XUUBEST(MAXNG),CORIND(MSTART)
	 CHARACTER*18 HS(7)
	 CHARACTER*14 H1(2)
         CHARACTER*33 NAMBEST
	 DATA (HS(J),J=1,7)/'NEAREST NEIGHBOUR ','FARTHEST NEIGHBOUR',
     &   'GROUP AVERAGE    ', 'MEDIAN            ','CENTROID          ',
     &   'FLEXIBLE SORTING  ','INCREMENTAL SUM SQR '/
         DATA (H1(J),J=1,2)/'UNSTANDARDIZED','STANDARDIZED  '/
         IF (FLAGS(11).EQ.0) WRITE(FYLENO,205)
         NH=0
205      FORMAT(/8X,'     Hierarchical-based initial groupings',/
     &          ,8X,'     ~~~~~~~~~~~~~~~~~~ ~~~~~~~ ~~~~~~~~~')
         OPEN(UNIT=28,FILE='hier.inp',STATUS='old',ERR=220)

210   CONTINUE
	 NH=NH+1
	 READ(28,*) ISU(NH),IS(NH)
	 IF (ISU(NH).EQ.-1) THEN
           CLOSE(28) 
           GOTO 230
         ENDIF
	 IF (IS(NH).EQ.6) THEN
	    READ(28,*) BETA(NH)
         ENDIF
      GOTO 210
 
220       CONTINUE
C       If input file hier.inp is not present these defaults are used
	  NH=7
	DATA (ISU(J),J=1,6)  /1,2,1,2,1,2/
	DATA (IS(J),J=1,6)  /3,3,2,2,7,7/

230      CONTINUE
	 NH=NH-1
         IF (FLAGS(11).EQ.0) THEN
          WRITE(FYLENO,*) '   ',NH,' Hierarchical Starts used'
         ENDIF
	 IF (NH.EQ.0) GOTO 260
         IF (FLAGS(3).EQ.1) THEN
  	   WRITE(*,*) 'Examining Hierarchical Starts...'
         ENDIF
c                   HIERARCHICAL STARTS
         DO 250 JJ=1,NH
c          Determine Hierarchical grouping
           CALL HIER(NIND,NATT,NG,X,IDT,ISU(JJ),IS(JJ),BETA(JJ),IFAULT)
	   IF (IFAULT.EQ.9) RETURN
c          Display Hierarchical grouping 
           IF (FLAGS(11).EQ.0) WRITE (FYLENO,235)
235        FORMAT (2X,'------------------------------------'
     &               ,'------------------------------------')
c           WRITE(FYLENO,*)'          ',JJ,' ',HS(IS(JJ))
c	   WRITE(FYLENO,*)'                ',H1(ISU(JJ))
            WRITE(FYLENO,*)'   ',JJ,' ',H1(ISU(JJ)),' ',HS(IS(JJ))
	    IF (IS(JJ).EQ.6) THEN
	       WRITE(FYLENO,*)'                BETA= ',BETA(JJ)
            ENDIF
          IF(FLAGS(11).EQ.0) THEN
             WRITE(FYLENO,245)(IDT(M),M=1,NIND)
245          FORMAT(2X,10I4)
             WRITE(FYLENO,*)
           ENDIF
           IF (FLAGS(3).EQ.1) THEN
            WRITE(*,*)'   ',JJ,' ',H1(ISU(JJ)),' ',HS(IS(JJ))
           ENDIF
C          Now we have the data and an IDT initial allocation
	   NS=NS+1
           CALL FIT(NIND,NATT,NG,NCOV,X,WL,XUU,XUUBEST,USA,
     &	    IDT,WBEST,XML,XMLBEST,NS,MAXNS,TOLS,W,CORIND,FFG,IER)
   	  WRITE(FYLENO,685) XML(NS)
685     FORMAT(2X,'Log Likelihood value from EM algorithm started',
     &        /2X,'from this grouping is',F15.3)
           IF (XML(NS).EQ.XMLBEST) THEN
             NAMBEST=H1(ISU(JJ))//' '//HS(IS(JJ))
           ENDIF
250    CONTINUE
260      RETURN
      END

      SUBROUTINE KMEANCON(NIND,NATT,NG,NCOV,X,WL,WBEST,XMLBEST,XML,
     &      NS,MAXNS,TOLS,NAMBEST,XUU,XUUBEST,USA,CORIND,FFG,IER)
C     This subroutine calls the Subroutine KMEANS and for the resulting
C     partition calls the subroutine FIT.
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
       INTEGER FLAGS(40),FYLENO,USA(MNIND)
       COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
       COMMON /STORE2/ FLAGS,FYLENO
       DIMENSION X(MNIND,MNATT),IDT(MNIND),WBEST(MNIND,MAXNG),
     &           XML(MSTART),TOLS(4),W(MNIND,MAXNG),CORIND(MSTART),
     &           WL(MNIND,MAXNG),XUU(MAXNG),XUUBEST(MAXNG)
         CHARACTER*33 NAMBEST
C       K-means Start
         IF (FLAGS(3).EQ.1) THEN
 	 WRITE(*,*) 'Examining K-means Starts'
           WRITE(*,*) 'K-means start Number',
     &     '  (Random seeds)'
           ENDIF
        DO 316 NKM=1,FLAGS(5)
        IF (FLAGS(11).EQ.0) THEN
        WRITE(FYLENO,305) NKM
305     FORMAT(2X,'           K-means Starting allocation',I3,/,
     &         2X,'           ~~~~~~~ ~~~~~~~~ ~~~~~~~~~~~')
        ENDIF
        IF (FLAGS(3).EQ.1) THEN
         IF (RANDTYPE.EQ.1.0) THEN
   	   WRITE(*,*) ' K-means Start',NKM,IX,IY,IZ
   	   WRITE(FYLENO,*) '   K-means Start',NKM,IX,IY,IZ
         ELSE
   	   WRITE(*,*) ' K-means Start',NKM,SEED
   	   WRITE(FYLENO,*) '  K-means Start',NKM,SEED
         ENDIF
        ENDIF
        CALL KMEANS(NIND,NATT,NG,X,IDT,EPSILON,MKMEAN,IER)
        IF (IER.EQ.-41) THEN 
        WRITE (FYLENO,*)'  WARNING : K-means did not converge current'
        WRITE (FYLENO,*)'            estimates will be used.'
         IER=0
        ENDIF
         IF(FLAGS(11).EQ.0) THEN
          WRITE(FYLENO,315)(IDT(MM),MM=1,NIND)
315       FORMAT(2X,10I4)
          WRITE(FYLENO,*)
         ENDIF
        NS=NS+1
        CALL FIT(NIND,NATT,NG,NCOV,X,WL,XUU,XUUBEST,USA,
     &	   IDT,WBEST,XML,XMLBEST,NS,MAXNS,TOLS,W,CORIND,FFG,IER)
   	  WRITE(FYLENO,685) XML(NS)
685     FORMAT(2X,'Log Likelihood value from EM algorithm started',
     &        /2X,'from this grouping is',F15.3)
        IF (XML(NS).EQ.XMLBEST) THEN
           I100=INT(FLOAT(NKM)/100.0)
           I10=INT(FLOAT(NKM-(I100*100))/10.0)
           I1=NKM-I100*100-I10*10
           NAMBEST='K-Means  '//CHAR(I100+48)
     &                            //CHAR(I10+48)
     &                            //CHAR(I1+48)
        ENDIF
316    CONTINUE
       RETURN
      END

       SUBROUTINE CHECK(NG,NIND,NATT,IDT,SWAPID)
C      This subroutine checks to see if the partition found has enough 
C      points in all of it's groups. If this is not so then the
C      subroutine swap is called and points are taken from other groups.

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
         DIMENSION IDT(MNIND),NO(MAXNG)
         DO 410 IF=1,NG
410        NO(IF)=0.0
         DO 420 IG=1,NIND
           IG1=IDT(IG)
           IF (IG1.GT.0) NO(IG1)=NO(IG1)+1
420      CONTINUE
         DO 430 IF1=1,NG
           IF(NO(IF1).LT.NATT+1) THEN
             CALL SWAP(IDT,NO,NIND,NATT,NG,IF1)
             SWAPID=IF1
           ENDIF
430     CONTINUE
        RETURN
       END


       SUBROUTINE SWAP(IDT,NO,NIND,NATT,NG,I)
C     This subroutine moves points from one group to another
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
         DIMENSION IDT(MNIND),NO(MAXNG)
        IDIF=NATT+2-NO(I)
        JJ=0
        DO 510 J=1,NG
          NOGIVE=NO(J)-NATT-2
          IF (NOGIVE.LT.0) GO TO 510
          IF(JJ.GE.IDIF)GO TO 510
          II=0
          DO 520 IG=1,NIND
            IF(II.GE.NOGIVE)GO TO 520
            IF(JJ.GE.IDIF)GO TO 510
            IF(IDT(IG).EQ.J) THEN
            IDT(IG)=I
            II=II+1
            JJ=JJ+1
            NO(J)=NO(J)-1
            NO(I)=NO(I)+1
          ENDIF
520      CONTINUE
510     CONTINUE
        RETURN
        END

      SUBROUTINE FIT(NIND,NATT,NG,NCOV,X,WL,XUU,XUUBEST,USA,
     &    IDT,WBEST,XML,XMLBEST,NS,MAXNS,TOLS,W,CORIND,FFG,IER)
C     This subroutine fits a mixture model to the sample contained in X
C     via the EM algorithm (subroutine LOOP) started from the initial
C     partition contained in IDT

       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
       INTEGER FLAGS(40),FYLENO,USA(MNIND)
       COMMON /STORE2/ FLAGS,FYLENO
       DIMENSION X(MNIND,MNATT),XVAR(MAXNG,MNATT,MNATT),
     &           T(MAXNG),DV(MAXNG),V(MAXNG,MNATT,MNATT),
     &           XMU(MAXNG,MNATT),XML(MSTART),IDT(MNIND),
     &           WBEST(MNIND,MAXNG),WL(MNIND),CORIND(MSTART),
     &           TOLS(4),W(MNIND,MAXNG),XUU(MAXNG),
     &           XUUBEST(MAXNG),TXUU(MAXNG)
           SWAPID=0
           CALL CHECK(NG,NIND,NATT,IDT,SWAPID)
	   IF (SWAPID.GT.0) THEN
	   WRITE(FYLENO,*)'  This grouping has too few points'
           WRITE(FYLENO,*)'  in group ',SWAPID,
     &     ' points will be re-allocated'
	   WRITE(FYLENO,*)'  from another group'
	    ENDIF
           IF ((SWAPID.GT.0).AND.(FLAGS(11).EQ.0)) THEN
             WRITE (FYLENO,*) '  New grouping is:'
	     WRITE(FYLENO,645)(IDT(MM),MM=1,NIND)
645          FORMAT(2X,10I3)
           ENDIF
      CALL ESTIMATES(NIND,NATT,NG,X,IDT,WL,NCOV,XMU,V,
     &	 XVAR,DV,T,USA,IER)
         IF (IER.NE.0) GOTO 660 
c	 ENDIF
           FLAGS(9)=2
        IF (FLAGS(7).EQ.4) THEN
          DO 10101 III=1,NIND
            DO 10100 KKK=1,NG
              W(III,KKK)=0
10100       CONTINUE
            IF (IDT(III).GT.0) THEN
              W(III,IDT(III))=1
            ENDIF
10101     CONTINUE

          CALL  TMOM(NIND,NATT,NG,X,XMU,XVAR,W,TXUU)
        ELSE
         DO 1999 II=1,NG
          TXUU(II)=XUU(II)
1999     CONTINUE
        ENDIF

	 CALL LOOP(NIND,NATT,NG,X,XMU,V,XVAR,DV,T,NCOV,IER,TXML,
     &          IDT,WL,W,TXUU,USA,TOLS)
	 IF (NS.GT.MAXNS) THEN
	   NS=MAXNS
	 ENDIF
	   XML(NS)=TXML
C         IF ((IER.EQ.2).AND.(XML(NS).EQ.0)) THEN
         IF (IER.EQ.2) THEN
           WRITE(FYLENO,*) '  No solution found'
           WRITE(FYLENO,*) '  This start will be ignored'
           WRITE(FYLENO,*) '  Log Like set to Zero'
           XML(NS)=0
	   CORIND(NS)=-1
           GOTO 700
         ENDIF
C         IF ((IER.EQ.5).AND.(XML(NS).EQ.0)) THEN
         IF (IER.EQ.5) THEN
           WRITE(FYLENO,*) '  No solution found'
           WRITE(FYLENO,*) '  This start will be ignored'
           WRITE(FYLENO,*) '  Log Like set to Zero'
           XML(NS)=0
	   CORIND(NS)=-1
           GOTO 700 
         ENDIF
   	 IF(IER-1)680,660,670
C        If the matrix has zero pivot set Log likelihood to 0 
C        This corresponds to convergence to a singular matrix
660      CONTINUE 
	 IF (NS.GT.MAXNS) THEN
	   NS=MAXNS
	 ENDIF
         XML(NS)=0.0
	 CORIND(NS)=-1
         WRITE(FYLENO,*) '   Problem with Inversion.'
         WRITE(FYLENO,*) '   CODE ',IER
         WRITE(FYLENO,*) '   Log likelihood set to Zero'
	 GOTO 680

670      IER=2
	 GOTO 700 

680      CONTINUE

c680	  WRITE(FYLENO,685) XML(NS)
c685     FORMAT(2X,'Log likelihood value for this grouping is',F15.3)
c         IF (XML(FLAGS(26)+NH+1).LT.XMLBEST) GOTO 583
         IF ((FFG.LT.0).AND.(CORIND(NS).GE.0)) THEN
          XMLBEST=XML(NS)
          FFG=1
         DO 1690 KK=1,NG
           XUUBEST(KK)=TXUU(KK)
           DO 1690 II=1,NIND
1690          WBEST(II,KK)=W(II,KK)

         ELSEIF ((XML(NS).GT.XMLBEST).AND.(CORIND(NS).GE.0)) THEN 
          XMLBEST=XML(NS)
          DO 690 KK=1,NG
            XUUBEST(KK)=TXUU(KK)
            DO 690 II=1,NIND
              WBEST(II,KK)=W(II,KK)
690       CONTINUE 
699      ENDIF 


700      CONTINUE
        RETURN
       END
C
C
C
       SUBROUTINE RESAMP(NIND,NATT,NG,NCOV,TX,TXMU,TXVAR,DVS,
     &    TT,TIDT,NREP,SEEDS,LLAM,TXML,TOLS,USA,
     &    SET,SEU,SEV,FSEED,XUU,IERS,IER)
C     Main subroutine for the program MMresamp written by D.Peel
C     Oct 1994- Oct1995 based on the program resamp.c by P.Adams 
C     Feb 1992

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     see head of driver file regarding this line
       INCLUDE 'EMMIX.max'
      EXTERNAL NMM,BSAMP,AUTOPARTITION
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM
      REAL LLAM
      INTEGER FLAGS(40),FYLENO
      INTEGER USA(MNIND),TIDT(MNIND)
      COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
      COMMON /STORE2/ FLAGS,FYLENO

      DIMENSION X(MNIND,MNATT),TX(MNIND,MNATT),WL(MNIND),
     &    XVAR(MAXNG,MNATT,MNATT),SEIER(MAXREP),
     &    T(MAXNG), XMU(MAXNG,MNATT),IDT(MNIND),LLAM(MAXREP),
     &    TXML(2,MAXREP),SEEDS(MAXREP,3),IERS(MAXREP,2),
     &    DVS(MAXREP,MAXNG),DV(MAXNG),TOLS(4),FSEED(3),XUU(MAXNG)
      DIMENSION TXVAR(MAXNG,MNATT,MNATT),TT(MAXNG),TXMU(MAXNG,MNATT),
     &    W(MNIND,MAXNG),XCC(MAXNG),TW(MNIND,MAXNG)
      DIMENSION SEEXPT(MAXNG),SEEXPU(MAXNG,MNATT),
     &    SEEXPV(MAXNG,MNATT,MNATT),
     &    SEMEAN(MAXREP,MAXNG,MNATT),SEMIX(MAXREP,MAXNG),
     &    SEVAR(MAXREP,MAXNG,MNATT,MNATT),SET(MAXNG),SEU(MAXNG,MNATT),
     &    SEV(MAXNG,MNATT,MNATT),V(MAXNG,MNATT,MNATT)
      DIMENSION DEN(MNIND,MAXNG+1),XLOGL(MITFIN),XMAH(MNIND,MAXNG)

      FYLENO=42
      IX=SEEDS(1,1)
      IY=SEEDS(1,2)
      IZ=SEEDS(1,3)

      IF ((FLAGS(18).GE.1).AND.(FLAGS(18).NE.4)) THEN
C       Obtain inverse and determinant of each estimated covariance
C       matrix
        CALL GDET(NATT,NG,TXVAR,V,DV,IER,NULL)
        IF (IER.NE.0) THEN
          WRITE (FYLENO,*) ' Problem inverting covariance matrix'
          WRITE (FYLENO,*) ' from current solution to calculate'
          WRITE (FYLENO,*) ' standard errors.'
          RETURN
        ENDIF
        IER=0
        DO 9 II=1,MNIND
          WL(II)=1
9       CONTINUE 
        IOUNT=1
        CALL ESTEP(NIND,NATT,NG,TX,TXMU,V,TT,DEN,WL,TW,XUU,USA,DV,
     &            XLOGL,IOUNT,XMAH,IER)
       ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Main Loop over replications
       DO 880 I=1,NREP
	 IERS(I,1)=0
         IER=0

C        Display progress to screen
         IF (FLAGS(10).EQ.0) THEN
          WRITE (*,*) '    Simulation',I,'/',NREP 
	  WRITE (*,804) SEEDS(I,1),SEEDS(I,2),SEEDS(I,3)
         ELSE
          WRITE (*,802) I,NREP,NG,NG+1
802       FORMAT(4X,'Simulation',I4,'/',I4,': ',I3,' vs ',I3,' groups')
	  WRITE (*,804) SEEDS(I,1),SEEDS(I,2),SEEDS(I,3)
804	  FORMAT(3X,'Seeds= ',F7.0,' ',F7.0,' ',F7.0)
         ENDIF
         WRITE(*,805)
805      FORMAT('  ------------------------------------------------')
         IF ((FLAGS(18).GT.0).AND.(FLAGS(18).NE.4)) THEN
         DO 806 II=1,NIND
           DO 806 JJ=1,NATT
             X(II,JJ)=TX(II,JJ)
806      CONTINUE
         ENDIF

C        Standard Errors
         IF (((FLAGS(18).NE.0).AND.(FLAGS(18).NE.4)).OR.
     &       (FLAGS(10).EQ.0)) THEN
         OPEN(UNIT=42,FILE='respSE.out',STATUS='UNKNOWN')
       	  OPEN(UNIT=43,FILE='SEsamp.out',STATUS='UNKNOWN')
	  FYLENO=43
          BT=0
          CALL BSAMP(NIND,NATT,NG,TXMU,TXVAR,TT,TX,X,WL,IDT,BT,IER)
	  CLOSE(43)
          
          IF (FLAGS(18).NE.0) THEN 
          DO 807 II=1,NIND
            DO 808 JJ=1,NG
C               W(II,JJ)=0.0
                W(II,JJ)=TW(II,JJ)
808         CONTINUE
C            W(II,TIDT(II))=1
807       CONTINUE
          ENDIF

          SEIER(I)=IER
          IF (IER.NE.0) GOTO 717 

C         Fit a Multivariate mixture model under Ho:No groups = NG to the
C         generated sample (X) using the true partition (IDT) from mvng
C         to start.
          DO 715 II=1,NIND
           IDT(II)=TIDT(II)
715       CONTINUE
          FYLENO=42
          
          IF (FLAGS(18).EQ.0) THEN
            DO 2816 KK=1,NG
            T(KK)=TT(KK)
            DO 2816 II=1,NATT
              XMU(KK,II)=TXMU(KK,II)
              DO 2816 JJ=1,NATT
                XVAR(KK,II,JJ)=TXVAR(KK,II,JJ)
2816         CONTINUE
  	    FLAGS(3)=3 
          ELSE
    	    FLAGS(3)=2 
            DO 2715 II=1,NIND
              DO 2715 JJ=1,NATT
                X(II,JJ)=TX(II,JJ)
2715        CONTINUE
          ENDIF

          CALL NMM(NIND,NATT,NG,NCOV,IDT,W,X,WL,TXML(1,I),
     &	           DV,V,TOLS,XMU,XVAR,T,XUU,USA,FSEED,IER)
          SEIER(I)=IER

C         For Standard Errors store parameter estimates
717         CONTINUE
            DO 740 II=1,NG
                SEMIX(I,II)=T(II)
              DO 740 JJ=1,NATT
                SEMEAN(I,II,JJ)=XMU(II,JJ)
              DO 740 LL=1,JJ
                SEVAR(I,II,JJ,LL)=XVAR(II,JJ,LL)
740         CONTINUE           

	  DO 730 KK=1,NG 
730          DVS(I,KK)=DV(KK) 
          CLOSE(42)
         ENDIF

C        Parametric BOOTSTRAPPING
         IF (FLAGS(10).GT.0) THEN
         OPEN(UNIT=42,FILE='respH0.out',STATUS='UNKNOWN')
C         Generate a random sample for the parametric bootstrap
          BT=1
       	  OPEN(UNIT=43,FILE='bsamp.out',STATUS='UNKNOWN')
	  FYLENO=43
          CALL BSAMP(NIND,NATT,NG,TXMU,TXVAR,TT,TX,X,WL,IDT,BT,IER)
          IF (IER.NE.0) THEN
            TXML(1,I)=0
            TXML(2,I)=0
            IERS(I,1)=4
            IERS(I,2)=IER
C           Report to file the problem for this bootstrap replicate
C           Problem in generation of sample
            IF (FLAGS(11).NE.2) THEN
              WRITE (42,815) I,IER
              WRITE (43,815) I,IER
815           FORMAT(I4,4X,
     &       'ERR: Problem in generation of sample (code',I2,')')
            ENDIF
            GOTO 840
          ENDIF
	  CLOSE(43)


C         Fit a Multivariate mixture model under Ho:No groups = NG to the
C         generated sample (X) using the true partition (IDT) from mvng
C         to start.
          DO 816 KK=1,NG
            T(KK)=TT(KK)
            DO 816 II=1,NATT
              XMU(KK,II)=TXMU(KK,II)
              DO 816 JJ=1,NATT
                XVAR(KK,II,JJ)=TXVAR(KK,II,JJ)
816       CONTINUE
          FYLENO=42
	  FLAGS(3)=3 
           IF ((FLAGS(17).EQ.0).OR.(FLAGS(18).EQ.4)) THEN
            FLAGS(21)=(-1)*FLAGS(21)
           ENDIF
          CALL NMM(NIND,NATT,NG,NCOV,IDT,W,X,WL,
     &	   TXML(1,I),DV,V,TOLS,XMU,XVAR,T,XUU,USA,FSEED,IER)
           
           IF ((FLAGS(17).EQ.0).OR.(FLAGS(18).EQ.4)) THEN
            FLAGS(21)=(-1)*FLAGS(21)
           ENDIF

C         If Standard Errors are needed and haven't been already stored 
C         store parameter estimates
          IF ((FLAGS(17).GT.0).AND.(FLAGS(18).EQ.0)) THEN
            SEIER(I)=IER
            DO 40 II=1,NG
                SEMIX(I,II)=T(II)
              DO 40 JJ=1,NATT
                SEMEAN(I,II,JJ)=XMU(II,JJ)
              DO 40 LL=1,JJ
                SEVAR(I,II,JJ,LL)=XVAR(II,JJ,LL)
40          CONTINUE           
          ENDIF

C         Detect if any errors have occurred, store error code and continue
          IF (IER.GT.0) THEN
	    TXML(1,I)=0
	    TXML(2,I)=0
	    IERS(I,1)=5
            IERS(I,2)=IER
C           Report to file the problem for this bootstrap replicate
C           Problem with fit under H0'
            WRITE (25,825) I,IER
825	    FORMAT(I4,4X,
     &       'ERR: Problem with fit under H0 (code',I2,')')
            GOTO 840
	  ENDIF
	  DO 830 KK=1,NG 
830          DVS(I,KK)=DV(KK) 
          CLOSE(42)

C         If Bootstrapping -2log Lambda then fit under H1
          IF (FLAGS(10).GT.0) THEN 
C           Find an initial partition via various clustering methods
            OPEN(UNIT=42,FILE='respH1.out',STATUS='UNKNOWN')
            FLAGS(21)=(-1)*FLAGS(21)
            CALL AUTOPARTITION(NIND,NATT,NG+1,NCOV,X,
     &                 W,WL,TOLS,XUU,USA,FSEED,IER)
C           Detect if any errors have occurred, store error code and continue
            IF (IER.GT.0) THEN
              TXML(2,I)=0
              IERS(I,1)=6
              IERS(I,2)=IER
C             Report to file the problem for this bootstrap replicate
C             Auto start unable to find start
              IF (FLAGS(11).NE.2) THEN
                WRITE (25,835) I,TXML(1,I),IER
              ENDIF
835           FORMAT(I4,G15.3,1X,
     &        'ERR:Auto start unable to find a start(code',I2,')')
	      GOTO 840
            ENDIF
            IF (FLAGS(11).NE.2) THEN
              CALL CAPART(NIND,NATT,NG+1,W,IDT,USA,XCC)
              WRITE(FYLENO,*) '  Start found'
              WRITE(FYLENO,836) (IDT(MM),MM=1,NIND)
            ENDIF
836         FORMAT(2X,10I4)
C           Fit a Multivariate mixture model under H1:No groups = NG+1 to the
C           generated sample (X) using various clustering methods.
            FYLENO=42
	    FLAGS(3)=2
            CALL NMM(NIND,NATT,NG+1,NCOV,IDT,W,X,WL,
     &      TXML(2,I),DV,V,TOLS,XMU,XVAR,T,XUU,USA,FSEED,IER)

            FLAGS(21)=(-1)*FLAGS(21)
C           Detect if any errors have occurred, store error code and continue
            IF (IER.GT.0) THEN
              TXML(2,I)=0
C             Report to file the problem for this bootstrap replicate
C             Problem with fit under H1
              IF (FLAGS(11).NE.2) THEN
                WRITE (25,845) I,TXML(1,I),IER
              ENDIF
845           FORMAT(I4,G15.3,4X,
     &        'ERR: Problem with fit under H1(code',I2,')')
              IERS(I,1)=7
	      IERS(I,2)=IER
              GOTO 840
            ENDIF
	    DO 850 KK=1,NG+1
              DVS(I,KK)=DV(KK) 
850         CONTINUE
		
            CLOSE(42)
C           Calculate -2log(lambda)
            LLAM(I)=(-2)*(TXML(1,I)-TXML(2,I))
            IF (FLAGS(11).NE.2) THEN
	      IF (LLAM(I).GT.0) THEN
                WRITE (25,855) I,TXML(1,I),TXML(2,I),LLAM(I)
855	        FORMAT(I4,G15.3,G15.3,G15.3)
              ELSE
	        WRITE (25,865) I,TXML(1,I),TXML(2,I),LLAM(I)
865	        FORMAT(I4,G15.3,G15.3,G15.3,' < Problem')
              ENDIF
            ENDIF		
           ENDIF

840         CONTINUE
            IF ((FLAGS(11).NE.2).AND.(FLAGS(10).GT.0)) THEN
              WRITE(25,875) SEEDS(I,1),SEEDS(I,2),SEEDS(I,3)		
              WRITE (25,*) '     DET= ',(DVS(I,KK),KK=1,NG)
            ENDIF
875         FORMAT (4X,'seeds= ',G15.5,G15.5,G15.5)
          ENDIF

          IF ((FLAGS(17).GT.0).AND.(FLAGS(18).NE.4)) THEN
          IF (SEIER(I).NE.0) THEN
            WRITE (26,*) '              SIMULATION ',I,
     &                 ' Problem code:',SEIER(I)
           ELSE
            WRITE (26,*) '              SIMULATION ',I
           ENDIF
          DO 877 K=1,NG
          WRITE (26,*) 'U grp',K,(SEMEAN(I,K,J),J=1,NATT)
          WRITE (26,*) 'Covariance',K
          DO 876 J=1,NATT
            WRITE (26,*) (SEVAR(I,K,J,L),L=1,J)
876       CONTINUE
877       CONTINUE
          IF (NG.GT.2) THEN
            WRITE (26,*) 'Mix proportions= ',(SEMIX(I,K),K=1,NG-1)
          ENDIF
         ENDIF

         SEEDS(I+1,1)=IX
         SEEDS(I+1,2)=IY
         SEEDS(I+1,3)=IZ
880   CONTINUE

      IF ((FLAGS(17).GT.0).AND.(FLAGS(18).NE.4)) THEN
         IF (NREP.GT.1) THEN        
            DO 270 K=1,NG
              SEEXPT(K)=0
              SET(K)=0
              DO 270 J=1,NATT
                SEU(K,J)=0
                SEEXPU(K,J)=0
                DO 270 I=1,J
                  SEV(K,J,L)=0
                  SEEXPV(K,J,I)=0
270         CONTINUE

            DO 300 I=1,NREP
                 DO 300 K=1,NG
             SEEXPT(K)=SEEXPT(K)+SEMIX(I,K)/REAL(NREP)
                   DO 290 J=1,NATT

             SEEXPU(K,J)=SEEXPU(K,J)+SEMEAN(I,K,J)/REAL(NREP)
                   DO 290 L=1,J
             SEEXPV(K,J,L)=SEEXPV(K,J,L)+SEVAR(I,K,J,L)/REAL(NREP)
290                CONTINUE
300      CONTINUE

          DO 305 I=1,NREP
            DO 310 K=1,NG
              SET(K)=SET(K)+(SEMIX(I,K)-SEEXPT(K))*
     &                       (SEMIX(I,K)-SEEXPT(K))
             DO 310 J=1,NATT
               SEU(K,J)=SEU(K,J)+(SEMEAN(I,K,J)-SEEXPU(K,J))*
     &                         (SEMEAN(I,K,J)-SEEXPU(K,J))
             DO 310 L=1,J
          SEV(K,J,L)=SEV(K,J,L)+(SEVAR(I,K,J,L)-SEEXPV(K,J,L))*
     &                         (SEVAR(I,K,J,L)-SEEXPV(K,J,L))
310        CONTINUE
305      CONTINUE


            DO 315 K=1,NG
             SET(K)=SQRT(SET(K)/REAL(NREP-1))
             DO 315 J=1,NATT
               SEU(K,J)=SQRT(SEU(K,J)/REAL(NREP-1))
               DO 315 L=1,J
                 SEV(K,J,L)=SQRT(SEV(K,J,L)/REAL(NREP-1))
315         CONTINUE
          ENDIF
          ENDIF
      RETURN
      END

      SUBROUTINE CALINFO(NIND,NATT,NG,X,XMU,V,W,T,NCOV,IDT,SET,SEU,SEV)
c      This subroutine displays all the relevant information from the
c      EM algorithm to the be used by the program MMresamp.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      COMMON /STORE3/ OFYLE,INFYLE,P1FYLE,P2FYLE,P3FYLE,OFYLE2
      DIMENSION T(MAXNG),XMU(MAXNG,MNATT),V(MAXNG,MNATT,MNATT),
     &  X(MNIND,MNATT),W(MNIND,MAXNG),IDT(MNIND),
     &  SET(MAXNG),SEU(MAXNG,MNATT),SEV(MAXNG,MNATT,MNATT),
     &  XINFO(MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2,
     &  MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2),
     &  VINFO(MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2,
     &  MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2),
     &  H(MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2),
     &  XPROD(MNATT)
      CHARACTER INFYLE*20, OFYLE*20
      CHARACTER P1FYLE*20,P2FYLE*20,P3FYLE*20
      CHARACTER OFYLE2*20
      ISIZE=NG-1+NG*NATT+NG*NATT*(NATT+1)/2
      DO 50 J=1,ISIZE
        DO 50 I=1,ISIZE
         XINFO(I,J)=0.0
50    CONTINUE

      DO 900 I=1,NIND
      IMARK=0
      DO 100 K=1,NG-1
         IMARK=IMARK+1
         H(IMARK)=W(I,K)/T(K)-W(I,NG)/T(NG)
100   CONTINUE   

      DO 200 K=1,NG
         DO 150 II=1,NATT
         XPROD(II)=0
         DO 140 JJ=1,NATT
          XPROD(II)=XPROD(II)+V(K,II,JJ)*(X(I,JJ)-XMU(K,JJ))
140      CONTINUE
         IMARK=IMARK+1
         H(IMARK)=W(I,K)*XPROD(II)
150      CONTINUE
200   CONTINUE   

      DO 300 K=1,NG
      DO 300 JJ=1,NATT
      DO 300 II=1,JJ
         XKRON=0
         IF (II.EQ.JJ) XKRON=1
         XTEMP1=0
         XTEMP2=0
      DO 250 III=1,NATT
          XTEMP1=XTEMP1+(X(I,III)-XMU(K,III))*V(K,III,II)
          XTEMP2=XTEMP2+(X(I,III)-XMU(K,III))*V(K,III,JJ)
250   CONTINUE
         IMARK=IMARK+1
         H(IMARK)=0.5*W(I,K)*(2-XKRON)*
     &            (-V(K,II,JJ)+(XTEMP1*XTEMP2))
300   CONTINUE
       DO 400 II=1,IMARK
       DO 400 JJ=1,IMARK
       XINFO(II,JJ)=XINFO(II,JJ)+H(II)*H(JJ)
400   CONTINUE
900   CONTINUE

      CALL XINFOINV(ISIZE,NG,XINFO,VINFO,IER,NULL)
      DO 902 K=1,NG
        SET(K)=0.0
        DO 902 I=1,NATT
          SEU(K,I)=0.0
          DO 902 J=1,NG
            SEV(K,I,J)=0.0
902   CONTINUE
      IMARK=0
      DO 920 K=1,NG-1
        IMARK=IMARK+1
        SET(K)=VINFO(IMARK,IMARK)
920   CONTINUE
      DO 930 K=1,NG
        DO 930 J=1,NATT
        IMARK=IMARK+1
        SEU(K,J)=VINFO(IMARK,IMARK)
930   CONTINUE
      DO 940 K=1,NG
        DO 940 J=1,NATT
        DO 940 I=1,J
        IMARK=IMARK+1
        SEV(K,J,I)=VINFO(IMARK,IMARK)
940   CONTINUE
      RETURN
      END

      SUBROUTINE XINFOINV(ISIZE,NG,XINFO,VINFO,IER,NULL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION
     &  XINFO(MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2,
     &  MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2),
     &  VINFO(MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2,
     &  MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2),
     &  TRIANG((MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2)*
     &  (MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2)),
     &  ATRIANG((MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2)*
     &  (MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2)),
     &  TEMP((MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2)*
     &  (MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2)),
     &  CMY((MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2)*
     &  (MAXNG-1+MAXNG*MNATT+MAXNG*MNATT*MNATT/2))
      NNULL=0
      DO 805 I=1,ISIZE*ISIZE
805      ATRIANG(I)=0
        IT=0
        TOL=0.0
        DO 810 J=1,ISIZE
          TOL=TOL+SQRT(XINFO(J,J))
          DO 810 I=1,J
            IT=IT+1
            TRIANG(IT)=XINFO(I,J)
810     CONTINUE
C        TOL=(TOL/ISIZE)*0.000001
        TOL=(TOL/ISIZE)*0.000000000000001
        CALL SYMINV (TRIANG,ISIZE,ATRIANG,TEMP,NULL,IER,RMAX,CMY,TOL)
        IF (IER.GT.0) GO TO 899
        IF (NULL.NE.0) THEN
          WRITE (FYLENO,815)  NULL
815       FORMAT (/2X,'Rank deficiency of information matrix ',
     &                ' is ',I3)
          NNULL=NULL+1
        ENDIF
        IT=0
        PROD=1.0
        DO 820 J=1,ISIZE
          JJ=J*(J+1)/2
          PROD=PROD*CMY(JJ)*CMY(JJ)
          DO 820 I=1,J
            IT=IT+1
            VINFO(I,J)=ATRIANG(IT)
            VINFO(J,I)=VINFO(I,J)
820     CONTINUE
        WRITE (FYLENO,*)'  Determinant of the information matrix= ',PROD
        IF (PROD.EQ.0) THEN
         WRITE (FYLENO,*)'  Determinant is equal to zero for
     & Information matrix'
          IER=5
          RETURN
        ENDIF
      NULL=NNULL
      RETURN
895   FORMAT (/2X,'Terminal error in SYMINV for Info matrix ',
     &        'as IFAULT is ',I3)
899    WRITE (FYLENO,895) IER
      WRITE (*,895)  IER
      RETURN
      END

C
C
C
       SUBROUTINE NMM(NIND,NATT,NG,NCOV,IDT,W,X,WL,
     &     TXML,DV,V,TOLS,XMU,XVAR,T,XUU,USA,FSEED,IER)
C     The purpose of these subroutines to fit a mixture model consisting of NG
C     multivariate normal components
C     This is done via the EM algorithm.
C
C     For information about how to run and use this program see the
C     file nmm.readme.
C
C     D.Peel Nov 1995 Based on the program kmmu by K.Basford
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,TEMPNCOV
      INTEGER USA(MNIND)
          COMMON /STORE2/ FLAGS,FYLENO
      EXTERNAL AUTOPARTITION
      EXTERNAL RANDNUM,DETERRANDOM
      DOUBLE PRECISION RANDNUM
      DIMENSION X(MNIND,MNATT),XVAR(MAXNG,MNATT,MNATT),T(MAXNG),
     &         DV(MAXNG),V(MAXNG,MNATT,MNATT),XMU(MAXNG,MNATT),
     &         IDT(MNIND),TOLS(4),FSEED(3),XUU(MAXNG),
     &         W(MNIND,MAXNG),WL(MNIND),U(MNIND,MNATT),
     &         XMAH(MNIND,MAXNG),XCC(MAXNG)
C
      IER=0  

C     Calculate total data size when data is weighted ie weighted 
C     likelihood and sum of weights does not equal one.
      WTOT=0.0
      DO 70 I=1,NIND
        WTOT=WTOT+WL(I)
70    CONTINUE

C  If only one group then fit is straightforward
      IF (NG.EQ.1) THEN
C       Set partition to all points in a single group
        DO 15 I=1,NIND
         IDT(I)=1
15      CONTINUE
C       Set Flag for final EM run and equal covariances
        FLAGS(9)=1
        TEMPNCOV=NCOV
	NCOV=1
C       Calculate parameter estimates 
        CALL ESTIMATES(NIND,NATT,NG,X,IDT,WL,NCOV,XMU,V,XVAR,
     &               DV,T,USA,IER)
        IF (IER.GT.0) GOTO 60
cc	CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
        IF (FLAGS(7).EQ.4) THEN
         DO 20101 III=1,NIND
          DO 20100 KKK=1,NG
           W(III,KKK)=0
20100    CONTINUE
         W(III,IDT(III))=1 
20101    CONTINUE
         CALL  TMOM(NIND,NATT,NG,X,XMU,XVAR,W,XUU)
         WRITE (FYLENO,*) '  Moment estimates of NU are '
         WRITE (FYLENO,*) '   ',(XUU(KKK),KKK=1,NG)
         ENDIF

C       Execute the EM algorithm
        CALL LOOP(NIND,NATT,NG,X,XMU,V,XVAR,DV,T,NCOV,IER,
     &            TXML,IDT,WL,W,XUU,USA,TOLS)
C       Go to end of routine to output section
       GOTO 23 
      ENDIF


C    BOOTSTRAP fit from given posterior probabilities
      IF (FLAGS(3).EQ.2) THEN
C      Correct posterior probabilities for classified data
       CALL PARCORR(NIND,NATT,NG,USA,W)
C      Calculate parameters from posterior probabilities
       TEMP=0
       CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,
     &                 DV,XMU,WTOT,T,W,XUU,XMAH,TEMP,U,IER)
C      Display parameter estimates to output file if required
       IF (FLAGS(11).EQ.0) THEN
        CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
        IF (IER.GT.0) GOTO 60
       ENDIF  

      ELSEIF (FLAGS(3).EQ.3) THEN
       IF (FLAGS(11).EQ.0) THEN
         CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
         IF (IER.GT.0) GOTO 60
       ENDIF
       FLAGS(3)=2
C       Obtain inverse of this matrix
	  CALL GDET(NATT,NG,XVAR,V,DV,IER,NULL)
	  IF (IER.GT.0) RETURN 
	    DO 992 K=2,NG
	      DV(K)=DV(1)
    	      DO 992 J=1,NATT
               DO 992 I=1,J
		V(K,I,J)=V(1,I,J)
		V(K,J,I)=V(K,I,J)
992          CONTINUE

C    User supplied partition of the data
      ELSEIF (FLAGS(4).EQ.1) THEN
C       Calculate parameter estimates from allocation 
        CALL ESTIMATES(NIND,NATT,NG,X,IDT,WL,NCOV,XMU,V,XVAR,DV,
     &                  T,USA,IER)
C       Display parameter estimates to output file if required
        IF (FLAGS(11).EQ.0) THEN
  	 CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
         IF (IER.GT.0) GOTO 60
        ENDIF

C    Use automatic starts
      ELSEIF (FLAGS(4).EQ.3) THEN 
C       Generate maximum likelihood solution from various starts
	CALL AUTOPARTITION(NIND,NATT,NG,NCOV,X,
     &	 W,WL,TOLS,XUU,USA,FSEED,IER)
	IF (IER.EQ.2) GOTO 60
C       Correct posterior probabilities for classified data
        CALL PARCORR(NIND,NATT,NG,USA,W)
C       Calculate parameters from posterior probabilities
        TEMP=0
        CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,
     &                 DV,XMU,WTOT,T,W,XUU,XMAH,TEMP,U,IER)

C       Display parameter estimates to output file if required
        IF (FLAGS(11).EQ.0) THEN
 	 CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
        ENDIF
        IF ((FLAGS(3).EQ.1).AND.(FLAGS(11).EQ.0)) THEN
          CALL CAPART(NIND,NATT,NG,W,IDT,USA,XCC)     
          WRITE (*,*) '  Best Partition'
          WRITE (*,315) (IDT(JJ),JJ=1,NIND)
        ENDIF
        IF (FLAGS(11).EQ.0) THEN
         CALL CAPART(NIND,NATT,NG,W,IDT,USA,XCC)     
         WRITE (FYLENO,*)  
         WRITE (FYLENO,*) '  Starting Partition Found'
         WRITE (FYLENO,315) (IDT(JJ),JJ=1,NIND)
315      FORMAT (2X,10I4)
         WRITE(FYLENO,*)
        ENDIF

C    User posterior probabilities (weights,fuzzy partition)  
      ELSEIF (FLAGS(4).EQ.4) THEN
C      Calculate parameters from posterior probabilities
       TEMP=0
       CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,
     &                 DV,XMU,WTOT,T,W,XUU,XMAH,TEMP,U,IER) 

      ENDIF 
C     END of options


C     Set flag to switch the EM algorithm on to Final iterations
C     conditions eg stopping conditions, tolerances
      FLAGS(9)=1
C     Set flag to switch off the Stochastic EM
      FLAGS(2)=0

       IF (IER.GT.0) THEN 
         WRITE (FYLENO,*)'  Unable to continue covariance singular'
         WRITE (FYLENO,*)'  too few points in one of the components'
         RETURN
       ENDIF
       IF (FLAGS(7).EQ.4) THEN
        DO 10101 III=1,NIND
         DO 10100 KKK=1,NG
          W(III,KKK)=0
10100    CONTINUE
         W(III,IDT(III))=1 
10101    CONTINUE
   
        CALL  TMOM(NIND,NATT,NG,X,XMU,XVAR,W,XUU)
        WRITE (FYLENO,*) '  Moment estimates of NU are '
        WRITE (FYLENO,*) '   ',(XUU(KKK),KKK=1,NG)
       ENDIF

C    Call MAIN EM Algorithm loop
      CALL LOOP(NIND,NATT,NG,X,XMU,V,XVAR,DV,T,NCOV,IER,
     &          TXML,IDT,WL,W,XUU,USA,TOLS)

23     CONTINUE
      CALL FINOUT(NIND,NATT,NCOV,NG,XUU,XMU,XVAR,T,IDT,XCC)
      IF (NG.EQ.1) NCOV=TEMPNCOV
      GOTO 99

60    WRITE (FYLENO,65)
65    FORMAT (/2X,'Terminal error in SYMINV from input data')

99    RETURN 
      END


      SUBROUTINE FINOUT(NIND,NATT,NCOV,NG,XUU,XMU,XVAR,T,IDT,XCC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      EXTERNAL AUTOPARTITION
      EXTERNAL RANDNUM,DETERRANDOM
      DOUBLE PRECISION RANDNUM
      DIMENSION XVAR(MAXNG,MNATT,MNATT),T(MAXNG),
     &         XMU(MAXNG,MNATT),IDT(MNIND),XUU(MAXNG),N(MAXNG),
     &         XCC(MAXNG)
C
C      FINAL OUTPUT Section

C    Output parameter estimates to file
      IF (FLAGS(7).GT.0) THEN
      WRITE (FYLENO,26)
26    FORMAT (/2X,'Estimated Nu for each component')
        WRITE (FYLENO,35) (XUU(K),K=1,NG)
      ENDIF
      IF (FLAGS(20).EQ.1) THEN
       DO 1030 K=1,NG
	WRITE (29,35) (XMU(K,J),J=1,NATT)
	IF (NCOV.EQ.1) THEN
	  DO 1040 J=1,NATT
            WRITE (29,45) (XVAR(1,I,J),I=1,J)
1040      CONTINUE
	ELSE
	  DO 1041 J=1,NATT
	    WRITE (29,45) (XVAR(K,I,J),I=1,J)
1041      CONTINUE
	ENDIF
1030   CONTINUE

1031   CONTINUE
       WRITE (29,45) (T(I),I=1,NG)
      ENDIF
      DO 1180 K=1,NG
1180    N(K)=0
      DO 1185 I=1,NIND
        K=IDT(I)
        IF (K.EQ.0) GO TO 1185
          N(K)=N(K)+1
1185  CONTINUE

      IF (NG.GT.1) THEN
      WRITE (FYLENO,1187)
1187  FORMAT (/2X,'Number assigned to each component')
      WRITE (FYLENO,1189) (N(K),K=1,NG)
1189   FORMAT (2X,10I6)
      WRITE (FYLENO,1191)
1191   FORMAT (/2X,'Estimate of mixing proportion for each component')
      WRITE (FYLENO,1193) (T(K),K=1,NG)
1193   FORMAT (2X,10F7.3)
C     Compute estimates of correct allocation rates
      CC=0.0
      DO 1195 K=1,NG
        XCC(K)=XCC(K)/(NIND*T(K))
        CC=CC+T(K)*XCC(K)
1195  CONTINUE
      WRITE (FYLENO,1197)
1197   FORMAT (/2X,'Estimates of correct allocation rates for ',
     &        'each component')
      WRITE (FYLENO,1193) (XCC(K),K=1,NG)
      WRITE (FYLENO,1198) CC
1198   FORMAT (/2X,'Estimate of overall correct allocation rate ',
     &        F7.3)
      ENDIF

      WRITE (FYLENO,25)
25    FORMAT (/2X,'Estimated mean (as a row vector) for each component')
      DO 30 K=1,NG
30      WRITE (FYLENO,35) (XMU(K,J),J=1,NATT)
35      FORMAT (2X,5G13.5)
C     Test if a common covariance matrix is specified (NCOV = 1)
      IF ((NCOV.EQ.1).OR.(NCOV.EQ.3)) THEN
       IF (NCOV.EQ.1) THEN
	WRITE (FYLENO,37)
37      FORMAT (/2X,'Estimated common covariance matrix ')
       ELSE
	WRITE (FYLENO,38)
38      FORMAT (/2X,'Estimated common diagonal covariance matrix ')
       ENDIF
	DO 40 J=1,NATT
40        WRITE (FYLENO,45) (XVAR(1,I,J),I=1,J)
45      FORMAT (5X,5G14.6)
      ELSEIF ((NCOV.EQ.2).OR.(NCOV.EQ.4)) THEN
        DO 50 K=1,NG
         IF (NCOV.EQ.2) THEN
          WRITE (FYLENO,46) K
46        FORMAT (/2X,'Estimated covariance matrix for component ',I2)
         ELSE
          WRITE (FYLENO,48) K
48    FORMAT(/2X,'Estimated diagonal covariance matrix for component '
     &,I2)
         ENDIF
          DO 50 J=1,NATT
            WRITE (FYLENO,47) (XVAR(K,I,J),I=1,J)
47      FORMAT (5X,5G14.6)
50      CONTINUE
      ENDIF
99    RETURN 
      END

      SUBROUTINE ESTEP(NIND,NATT,NG,X,XMU,V,T,DEN,WL,W,XUU,USA,DV,
     &                 XLOGL,IOUNT,XMAH,IER)
C     This Subroutine implements the E-step of the EM algorithm 
      implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      INTEGER USA(MNIND)
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION X(MNIND,MNATT),XMU(MAXNG,MNATT),DV(MAXNG),
     &          T(MAXNG),V(MAXNG,MNATT,MNATT),
     &          WL(MNIND),W(MNIND,MAXNG),XUU(MAXNG),
     &          DEN(MNIND,MAXNG+1),XLOGL(MITER),AL(MAXNG),
     &          XMAH(MNIND,MAXNG)

      DO 901 K=1,NG
	IF (DV(K).EQ.0) THEN
	    IER=5
	    RETURN
        ENDIF
901   CONTINUE
       PI=3.141592653
       XLOGL(IOUNT)=0.0
C      Calculate multivariate density of each point for every group
       DO 950 JJ=1,NIND
        GUM=0.0
        DO 920 K=1,NG
          AL(K)=0.0
          XMAH(JJ,K)=0.0
          DO 910 I=1,NATT
            DO 910 J=1,NATT
              XMAH(JJ,K)=XMAH(JJ,K)+(X(JJ,I)-XMU(K,I))*V(K,I,J)*
     &              (X(JJ,J)-XMU(K,J))
              AL(K)=AL(K)+(X(JJ,I)-XMU(K,I))*V(K,I,J)*
     &	            (X(JJ,J)-XMU(K,J))
910       CONTINUE

          IF (FLAGS(7).GT.0) THEN
c     The NU is limited to stop floating point overflow
c     to increase this limit modify the 'EMMIX.max' file 
      IF (XUU(K).GT.NUMAX) THEN
	XUU(K)=NUMAX
        IER=-53
      ENDIF
c      IF (XUU(K).LT.300) THEN
      AL(K)=EXP(GAMMLN(XUU(K)/2.0+FLOAT(NATT)/2.0))
      AL(K)=AL(K)/(SQRT(DV(K)))
      AL(K)=AL(K)/(3.14*XUU(K))**(FLOAT(NATT)/2.0)
      AL(K)=AL(K)/EXP(GAMMLN(XUU(K)/2.0))
      AL(K)=AL(K)/(1.0+XMAH(JJ,K)/XUU(K))**(1/2.0*(XUU(K)+FLOAT(NATT)))
c      ELSE
c	XUU(K)=300
c        IER=-53
c        RETURN
c      ENDIF
          ELSE

C         Check as if AL(K) too large under flow may occur 
c          IF (AL(K).GT.75.0) THEN 
c           IF (AL(K).GT.175.0) THEN 
          IF (AL(K).GT.DENMAX) THEN 
            AL(K)=0.0
          ELSE
C           Calculate component density
            AL(K)=(-0.5)*AL(K)
            AL(K)=EXP(AL(K))/(SQRT(DV(K))*(2.*PI)**(NATT/2.0))
          ENDIF
         ENDIF
C        User (Pierre Prado) requested option
          IF (FLAGS(6).EQ.1) THEN
            DEN(JJ,K)=AL(K) 
          ENDIF
C         Check to stop underflow
          IF ((T(K).GT.XLOWEM).AND.(AL(K).GT.XLOWEM)) THEN 
C           Calculate mixture density            
            GUM=GUM+T(K)*AL(K)
          ELSE
          ENDIF
920     CONTINUE
C       Store mixture density
        DEN(JJ,NG+1)=GUM

C       Compute current estimates of posterior probabilities of group
C       membership (W)

C       Check if already classified then set to classified group 
        IF ((FLAGS(12).GT.0).AND.(USA(JJ).GT.-1)) THEN
          DO 59 J=1,NG
            W(JJ,J)=0
59        CONTINUE
          W(JJ,USA(JJ))=1
          IF ((AL(USA(JJ)).GT.XLOWEM).AND.(T(USA(JJ)).GT.XLOWEM)) THEN
           XLOGL(IOUNT)=XLOGL(IOUNT)+LOG(AL(USA(JJ))*T(USA(JJ)))*WL(JJ)
          ENDIF
        ELSE
        IF (GUM.EQ.0.0) THEN
          DO 930 K=1,NG
930         W(JJ,K)=0.0
	    IER=-111
        ELSE
          DO 940 K=1,NG
C           Check to catch numerical underflow
            IF (T(K).LT.XLOWEM.OR.AL(K).LT.XLOWEM) THEN 
              W(JJ,K)=0.0
            ELSE
C             Calculate posterior probabilities
              W(JJ,K)=T(K)*AL(K)/GUM
              W(JJ,K)=W(JJ,K)*WL(JJ)
            ENDIF
940       CONTINUE
C         Calculate Log-Likelihood contribution from point JJ
         
          XLOGL(IOUNT)=XLOGL(IOUNT)+LOG(GUM)*WL(JJ)
        ENDIF
       ENDIF
950   CONTINUE
      RETURN
      END 

      SUBROUTINE MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,
     &                 DV,XMU,WTOT,T,W,XUU,XMAH,TEMP,U,IER)
C     This Subroutine implements the M-step of the EM algorithm 
      implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION X(MNIND,MNATT),XMU(MAXNG,MNATT),
     &          XVAR(MAXNG,MNATT,MNATT),U(MNIND,MAXNG),
     &          T(MAXNG),V(MAXNG,MNATT,MNATT),DV(MAXNG),
     &          XV(MNATT,MNATT),WSUM(MAXNG),W(MNIND,MAXNG),
     &          WUSUM(MAXNG),XMAH(MNIND,MAXNG),XUU(MAXNG)
      IF ((FLAGS(7).GT.0).AND.(TEMP.EQ.1.0)) THEN
        DO 1370 K=1,NG
          WUSUM(K)=0
          WSUM(K)=0
          DO 1369 J=1,NIND
            U(J,K)=(XUU(K)+NATT)/(XUU(K)+XMAH(J,K))
            WUSUM(K)=WUSUM(K)+W(J,K)*U(J,K)
            WSUM(K)=WSUM(K)+W(J,K)
1369    CONTINUE
1370    CONTINUE
        IF (FLAGS(7).GT.1) THEN
          CALL TFREE(NIND,NATT,NG,XUU,U,W,1,IER)
        ENDIF
      ELSE
C     Compute new estimate of mixing proportion (T) for each group
      DO 960 K=1,NG
C       Calculate sum of posterior probabilities in each group
        WSUM(K)=0.0
        WUSUM(K)=0
        DO 960 JJ=1,NIND
          U(JJ,K)=1
          WSUM(K)=WSUM(K)+W(JJ,K)
          WUSUM(K)=WSUM(K)
960     CONTINUE
      ENDIF
      DO 970 K=1,NG
         IF (WUSUM(K).EQ.0) THEN
          IER=22
        WRITE(FYLENO,*)'  Problem: No points allocated to component ',K
        WRITE(FYLENO,*)'           during EM iteration '
        RETURN
         ENDIF
C       Calculate mixing proportion for group K
        T(K)=WSUM(K)/WTOT
970   CONTINUE
C     Compute new estimates of group means and covariance matrices
c       Note: matrix V is not inverse here but used to as temp storage
c             for un-inverted matrix
         CALL LCAL(NIND,NATT,NG,X,W,XMU,WSUM,WUSUM,XVAR,U)
C     Test if a common covariance matrix is specified (NCOV = 1)
C     or common covariance matrix with diag matrices (NCOV = 3)
      IF ((NCOV.EQ.1).OR.(NCOV.EQ.3)) THEN
C       Compute new estimate of common covariance matrix
        DO 990 J=1,NATT
          DO 990 I=1,J
            XV(I,J)=0.0
            DO 980 K=1,NG
980           XV(I,J)=XV(I,J)+XVAR(K,I,J)*WSUM(K)
              XVAR(1,I,J)=XV(I,J)/WTOT
              XVAR(1,J,I)=XVAR(1,I,J)
990     CONTINUE
        IF (NCOV.EQ.3) THEN
         DO 989 I=1,NATT
          DO 989 J=1,NATT
            IF (I.NE.J) XVAR(1,I,J)=0
989      CONTINUE
        ENDIF
        DO 991 K=2,NG
        DO 991 J=1,NATT
          DO 991 I=1,J
            XVAR(K,I,J)=XVAR(1,I,J)
            XVAR(K,J,I)=XVAR(K,I,J)
991     CONTINUE
       
C       Obtain inverse of this matrix
        CALL GDET(NATT,1,XVAR,V,DV,IER,NULL)
        IF (IER.GT.0) GO TO 1030
        DO 992 K=2,NG
          DV(K)=DV(1)
        DO 992 J=1,NATT
          DO 992 I=1,J
            V(K,I,J)=V(1,I,J)
            V(K,J,I)=V(K,I,J)
992     CONTINUE
      ELSEIF ((NCOV.EQ.2).OR.(NCOV.EQ.4)) THEN
         IF (NCOV.EQ.4) THEN
         DO 993 K=1,NG
         DO 993 I=1,NATT
          DO 993 J=1,NATT
            IF (I.NE.J) XVAR(K,I,J)=0
993      CONTINUE
         ENDIF
C       Obtain inverse and determinant of each estimated covariance
C       matrix
        CALL GDET(NATT,NG,XVAR,V,DV,IER,NULL)
        IF (IER.GT.0) GO TO 1030
      ENDIF
      IF (NULL.NE.0) THEN
        WRITE (FYLENO,1025)
1025    FORMAT (/2X,'Note: Because of possible problems with SYMINV ',
     &        'the program will stop ',/8X,'and print results from ',
     &        'current cycle estimates.')
        IER=4
        RETURN 
      ENDIF

      RETURN

1030      WRITE (FYLENO,1035)
1035  FORMAT (/2X,'Note: Because of problems with SYMINV this program',
     &       ' will ',/8X,'print results from current estimates.',/8X,
     &       'New means but old covariance matrices are printed.')
      RETURN
      END

      SUBROUTINE PARCORR(NIND,NATT,NG,USA,W)
C      This subroutine corrects the posterior probabilities
C      for classified data by reassigning the classified points 
C      to their prior groups
       implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX.max'
       INTEGER FLAGS(40),FYLENO,USA(MNIND)
       COMMON /STORE2/ FLAGS,FYLENO
       DIMENSION W(MNIND,MAXNG)
           DO 60 JJ=1,NIND
            IF ((FLAGS(12).GT.0).AND.(USA(JJ).GT.-1)) THEN
             DO 59 J=1,NG
              W(JJ,J)=0
59           CONTINUE
             W(JJ,USA(JJ))=1
            ENDIF
60          CONTINUE

         RETURN
        END


      SUBROUTINE ESTIMATES(NIND,NATT,NG,X,IDT,WL,NCOV,
     &                     XMU,V,XVAR,DV,T,USA,IER)
C     This Subroutine handles the estimation of the parameters
C     of a multivariate normal mixture model given a sample X and an
C     initial partition (group allocation) IDT.
      implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,USA(MNIND)
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION X(MNIND,MNATT),XMU(MAXNG,MNATT),U(MNIND,MAXNG),
     &          XVAR(MAXNG,MNATT,MNATT),IDT(MNIND),WUSUM(MAXNG),
     &          T(MAXNG),V(MAXNG,MNATT,MNATT),DV(MAXNG),
     &          XV(MNATT,MNATT),WL(MNIND),WSUM(MAXNG),W(MNIND,MAXNG)
C     Compute estimates of mean, covariance matrix and mixing
C     proportion for each group
      ICOUNT=0
      DO 50 K=1,NG
        DO 50 I=1,NIND
         U(I,K)=1
50    CONTINUE
      IF ((FLAGS(8).NE.4).AND.(FLAGS(12).GT.0)) THEN  
        DO 60 II=1,NIND
          IF (USA(II).GT.-1) THEN
            IDT(II)=USA(II)
            ICOUNT=ICOUNT+1
          ENDIF  
60      CONTINUE
      ELSEIF (FLAGS(12).GT.0) THEN  
        DO 61 II=1,NIND
          IF (USA(II).GT.-1) THEN
            ICOUNT=ICOUNT+1
            IDT(ICOUNT)=USA(II)
          ENDIF  
61      CONTINUE
	NTNIND=NIND
	NIND=ICOUNT
      ENDIF
      WTOT=0.0
      DO 70 I=1,NIND
        WTOT=WTOT+WL(I)
70    CONTINUE
      DO 80 K=1,NG
        WSUM(K)=0
80    CONTINUE
      DO 95 I=1,NIND
        DO 90 K=1,NG 
          W(I,K)=0
90      CONTINUE
        KK=IDT(I)
        IF (KK.NE.0) THEN
          W(I,KK)=WL(I)
          WSUM(KK)=WSUM(KK)+WL(I)
        ENDIF
95    CONTINUE
      DO 99 K=1,NG
         WUSUM(K)=WSUM(K)
         IF (WSUM(K).EQ.0) THEN
           IER=21
           RETURN
         ENDIF
99       T(K)=WSUM(K)/WTOT
       CALL LCAL(NIND,NATT,NG,X,W,XMU,WSUM,WUSUM,XVAR,U)
C     Test if a common covariance matrix is specified (NCOV = 1)
C     or common covariance matrix with diag matrices (NCOV = 3)
      IF ((NCOV.EQ.1).OR.(NCOV.EQ.3)) THEN
C       Compute pooled estimate of common covariance matrix
        DO 110 J=1,NATT
          DO 110 I=1,J
          XV(I,J)=0.0
          DO 100 K=1,NG
100         XV(I,J)=XV(I,J)+WSUM(K)*XVAR(K,I,J)
          XVAR(1,I,J)=XV(I,J)/(WTOT-FLOAT(NG))
          XVAR(1,J,I)=XVAR(1,I,J)
110     CONTINUE
        IF (NCOV.EQ.3) THEN
         DO 111 I=1,NATT
          DO 111 J=1,NATT
            IF (I.NE.J) XVAR(1,I,J)=0
111      CONTINUE
        ENDIF
        DO 112 K=2,NG
        DO 112 J=1,NATT
          DO 112 I=1,J
            XVAR(K,I,J)=XVAR(1,I,J)
            XVAR(K,J,I)=XVAR(K,I,J)
112     CONTINUE
C       Obtain inverse of this matrix
        CALL GDET(NATT,1,XVAR,V,DV,IER,NULL)
        IF (IER.GT.0) RETURN
        DO 120 K=2,NG
          DV(K)=DV(1)
          DO 120 J=1,NATT
            DO 120 I=1,NATT
              XVAR(K,I,J)=XVAR(1,I,J)
              V(K,I,J)=V(1,I,J)
120     CONTINUE
      ELSE
       IF (NCOV.EQ.4) THEN
        DO 121 K=1,NG
          DO 121 I=1,NATT
            DO 121 J=1,NATT
            IF (I.NE.J) XVAR(K,I,J)=0
121      CONTINUE
       ENDIF 
C       Obtain inverse and determinant of each estimated covariance
C       matrix
        CALL GDET(NATT,NG,XVAR,V,DV,IER,NULL)
      ENDIF
       IF ((FLAGS(8).EQ.4).AND.(FLAGS(12).GT.0)) THEN
	NIND=NTNIND
      ENDIF
      RETURN
      END
      
      SUBROUTINE OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
C     This subroutine writes the estimates of the parameters to the outfile
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
       INCLUDE 'EMMIX.max'
      DIMENSION XMU(MAXNG,MNATT),XVAR(MAXNG,MNATT,MNATT),
     &       T(MAXNG),V(MAXNG,MNATT,MNATT),DV(MAXNG),XUU(MAXNG)
      DO 210 K=1,NG
      WRITE (FYLENO,205) K
205   FORMAT (/2X,'Estimated mean (as a row vector) for component ',I2)
210     WRITE (FYLENO,215) (XMU(K,J),J=1,NATT)
215   FORMAT (2X,5G14.6)
      IF ((NCOV.EQ.2).OR.(NCOV.EQ.4)) THEN 
        DO 220 K=1,NG
         IF (NCOV.EQ.2) THEN
          WRITE (FYLENO,217) K
217       FORMAT (/2X,'Estimated covariance matrix for component ',I4)
         ELSE
          WRITE (FYLENO,231) K
231   FORMAT(/2X,'Estimated diagonal covariance matrix for component '
     & ,I4)
         ENDIF
          DO 220 J=1,NATT
            WRITE (FYLENO,218) (XVAR(K,I,J),I=1,J)
218         FORMAT (2X,5G14.6)
220     CONTINUE
      ELSE 
        IF (NCOV.EQ.1) THEN
         WRITE (FYLENO,225)
225      FORMAT (/2X,'Estimated common component-covariance matrix ')
        ELSE
         WRITE (FYLENO,226)
226      FORMAT(/2X,'Estimated common diagonal component-covariance',
     &   ' matrix')
        ENDIF
        DO 230 J=1,NATT
230       WRITE (FYLENO,218) (XVAR(1,I,J),I=1,J)
      ENDIF

      WRITE (FYLENO,235)
235   FORMAT (/2X,'Mixing proportion from each component')
      WRITE (FYLENO,237) (T(K),K=1,NG)
237   FORMAT (5X,10F7.3)

      IF (FLAGS(7).GT.0) THEN
        WRITE (FYLENO,*) 'Initial estimates of Nu'
        WRITE (FYLENO,*) (XUU(K),K=1,NG) 
      ENDIF
      RETURN
      END


      SUBROUTINE GDET(NATT,NG,XVAR,V,DV,IER,NULL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine reads all covariance matrices, then calls
C     SYMINV, which inverts a matrix and calculates its determinant,
C     for each covariance matrix in turn.
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION XVAR(MAXNG,MNATT,MNATT),V(MAXNG,MNATT,MNATT),
     &          DV(MAXNG),
     &          TRIANG(MNATT*MNATT),ATRIANG(MNATT*MNATT),
     &          TEMP(MNATT*MNATT),CMY(MNATT*MNATT)
      NNULL=0
      DO 805 I=1,NATT*NATT
805      ATRIANG(I)=0
      DO 830 K=1,NG
        IT=0
        TOL=0.0
        DO 810 J=1,NATT
          TOL=TOL+SQRT(XVAR(K,J,J))
          DO 810 I=1,J
            IT=IT+1
            TRIANG(IT)=XVAR(K,I,J)
810     CONTINUE
        TOL=(TOL/NATT)*0.000001
        CALL SYMINV (TRIANG,NATT,ATRIANG,TEMP,NULL,IER,RMAX,CMY,TOL)
        IF (IER.GT.0) GO TO 899 
        IF (NULL.NE.0) THEN
          WRITE (FYLENO,815)  K, NULL
815       FORMAT (/2X,'Rank deficiency of covariance matrix ',I3,
     &                ' is ',I3)
          
cc          WRITE (*,815) K, NULL
          NNULL=NULL+1
        ENDIF
        IT=0
        PROD=1.0
        DO 820 J=1,NATT
          JJ=J*(J+1)/2
          PROD=PROD*CMY(JJ)*CMY(JJ)
          DO 820 I=1,J
            IT=IT+1
            V(K,I,J)=ATRIANG(IT)
            V(K,J,I)=V(K,I,J)
820     CONTINUE
        DV(K)=PROD
        IF (DV(K).EQ.0) THEN
         WRITE (FYLENO,*)'  Determinant is equal to zero for grp ',K
          IER=5
          RETURN
        ENDIF
830   CONTINUE
      NULL=NNULL
      RETURN

895   FORMAT (/2X,'Terminal error in SYMINV for matrix ',I3,
     &        'as IFAULT is ',I3)
899    WRITE (FYLENO,895) K, IER
      WRITE (*,895)  K, IER
      RETURN
      END


      SUBROUTINE LOOP(NIND,NATT,NG,X,XMU,V,XVAR,DV,
     &                T,NCOV,IER,TXML,IDT,WL,W,XUU,USA,TOLS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     This subroutine uses the EM algorithm from a specified starting
C     value to find a solution of the likelihood equation.
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,USA(MNIND)
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION X(MNIND,MNATT),V(MAXNG,MNATT,MNATT),
     &          DV(MAXNG),T(MAXNG),W(MNIND,MAXNG),
     &          XUU(MAXNG),XMAH(MNIND,MAXNG),
     &          XVAR(MAXNG,MNATT,MNATT),IDT(MNIND),TOLS(4),
     &          XLOGL(MITER),WL(MNIND),XCC(MAXNG),
     &          XMU(MAXNG,MNATT),DEN(MNIND,MAXNG+1),
     &          U(MNIND,MAXNG),XLA(MITER)

c     initialise some parameters
       WTOT=0.0
       XLANEW=-10000000
       IOUNT=1
c     calculate the data size when data is weighted
       DO 979 I=1,NIND
         WTOT=WTOT+WL(I)
         DO 979 K=1,NG
979    CONTINUE

c        IF (FLAGS(7).EQ.4) THEN
c        DO 10101 III=1,NIND
c         DO 10100 KKK=1,NG
c          W(III,KKK)=0
c10100    CONTINUE
c         W(III,IDT(III))=1 
c10101    CONTINUE
c        CALL  TMOM(NIND,NATT,NG,X,XMU,XVAR,W,XUU)
c        WRITE (FYLENO,*) '  Moment estimates of NU are '
c        WRITE (FYLENO,*) '   ',(XUU(KKK),KKK=1,NG)
c        ENDIF

c     Set up stopping criteria either to default or user defined
        IF (FLAGS(9).EQ.1) THEN
          TOLERENCE=TOLS(3)
	  MAXIT=TOLS(4)
        ELSE
          TOLERENCE=TOLS(1)
	  MAXIT=TOLS(2)
        ENDIF

c     initialise weights to zero
       DO 901 I=1,NIND
 	 DO 901 J=1,NG
  	  W(I,J)=0.0
901    CONTINUE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       MAIN ITERATIVE LOOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
900   CONTINUE

C     Display to screen determinants if required 
       IF ((FLAGS(9).EQ.1).AND.(NG.NE.1)
     &     .AND.(FLAGS(11).EQ.0).AND.(FLAGS(3).EQ.1)) THEN
         WRITE (FYLENO,905)
         IF (FLAGS(8).EQ.2) WRITE (*,905)
905      FORMAT (/2X,
     &   'Determinants of component covariance matrices')
         WRITE (FYLENO,*) (DV(K),K=1,NG)
         IF (FLAGS(8).EQ.2) WRITE (*,*) (DV(K),K=1,NG)
       ENDIF
       CALL ESTEP(NIND,NATT,NG,X,XMU,V,T,DEN,WL,W,XUU,USA,DV,
     &            XLOGL,IOUNT,XMAH,IER)
c       IF (IER.GT.0) GOTO 9999 
       IF (IER.GT.0) GOTO 1099 

       IF (FLAGS(2).EQ.1) THEN
         CALL SEM(NIND,NATT,NG,W)
       ENDIF
C       CALL CAPART(NIND,NATT,NG,W,IDT,USA,XCC)     
C       WRITE (*,11111) (IDT(III),III=1,NIND)
C11111  FORMAT(10I3)
C     Display to screen Log-likelihood if required
       IF ((FLAGS(9).EQ.1).AND.(NG.NE.1)
     &     .AND.(FLAGS(11).EQ.0)) THEN
         WRITE (FYLENO,955) IOUNT-1, XLOGL(IOUNT)
         IF (FLAGS(21).LT.0) THEN
            WRITE (FYLENO,956) XLA(IOUNT)
         ENDIF
         IF ((FLAGS(8).EQ.2).AND.(FLAGS(3).EQ.1)) THEN
           WRITE (*,*) '  After iteration ',IOUNT-1,
     &     ' the likelihood = ',XLOGL(IOUNT)
955   FORMAT (2X,'After iteration ',I3,' the log likelihood =  ',F15.3)
956   FORMAT (2X,'  Aitkens accelerated estimate =  ',F15.3)
	   IF (IER.EQ.-53) THEN
	  WRITE(FYLENO,*)'  Warning: Estimated Nu exceeds 300'
	  WRITE(FYLENO,*)'           (Nu set to equal 300)'
	   ENDIF
	   IF (IER.EQ.-111) THEN
	  WRITE(FYLENO,*)'  Warning: Some points have zero Likelihood'
	  WRITE(FYLENO,*) '         (will denote with 0 in grouping)'
	  WRITE(*,*) 'Warning : Some points have zero Likelihood'
	   ENDIF
         ENDIF
       ENDIF


C     Test for exit from loop
       IF (IOUNT.GE.MAXIT) THEN 
         WRITE (FYLENO,115) MAXIT
115      FORMAT (/2X,'Note: This sample did not converge in ',I3,
     &     ' iterations.',/8X,'However the program will continue ',
     &     'to print results ',/8X,'obtained from the last cycle ',
     &     'estimates.')
         GO TO 1099
       ENDIF
C      Standard stopping criteria
       IF (IOUNT.GT.10) THEN 
        LAST=IOUNT-10
C       Aitkin's acceleration to be used when specified and doing
C       a bootstrap fit
ccc        IF (FLAGS(21).GE.0) THEN
         ALIM=TOLERENCE*XLOGL(LAST)
         DIFF=XLOGL(IOUNT)-XLOGL(LAST)
         IF (ABS(DIFF).LE.ABS(ALIM)) THEN
           XLA(IOUNT)=0
           GO TO 1099 
         ENDIF


         IF (FLAGS(21).LT.0) THEN
ccc        ELSE
          XLAOLD=XLA(LAST)
          XNUM=XLOGL(IOUNT)-XLOGL(IOUNT-1)
          DEM=XLOGL(IOUNT-1)-XLOGL(IOUNT-2)
         IF (DEM.LT.1E-35) THEN
          XLANEW=0
         ELSE
        C=XNUM/DEM
 
         ENDIF
         IF ((C.LT.(.99)).OR.(C.GT.(1.01))) THEN
          XLANEW=XLOGL(IOUNT-1)+XNUM*1/(1-C)
         ELSE
          XLANEW=0
         ENDIF
 
          XLA(IOUNT)=XLANEW
 
         IF (XLA(IOUNT).NE.0) THEN
          ALIM=TOLERENCE*XLAOLD
          DIFF=XLA(IOUNT)-XLAOLD
         ENDIF
      IF ((ABS(DIFF).LE.ABS(ALIM)).AND.(XLA(IOUNT).GE.XLOGL(IOUNT)))
     &      THEN
          XLOGL(IOUNT)=XLA(IOUNT)        
          GOTO 1099 
         ENDIF
       ENDIF
      ENDIF

      TEMP=1.0
c       IF (IOUNT.EQ.1) THEN
c       DO 1002 III=1,NIND
c       DO 1002 KKK=1,NG
c         U(III,KKK)=1
c1002   CONTINUE
c       CALL TFREE(NIND,NATT,NG,XUU,U,W,0,IER)
c       ENDIF
      
      CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,DV,
     &           XMU,WTOT,T,W,XUU,XMAH,TEMP,U,IER)
cccc      IF (IER.GT.0) GOTO 9999 
      IF (IER.GT.0) GOTO 1099 
      IOUNT=IOUNT+1
      GOTO 900
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     END OF MAIN ITERATIVE LOOP
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




C    Final work after convergence
1099  CONTINUE
      IF (IER.EQ.-111) THEN
	WRITE(*,*) 'WARNING : Some points have zero Likelihood'
        WRITE(*,*) '         (will denote with 0 in grouping)'
        IER=0
      ENDIF
      IF ((FLAGS(21).LT.0).AND.(XLA(IOUNT).NE.0).AND.
     &  (XLA(IOUNT).GE.XLOGL(IOUNT))) THEN
       TXML=XLA(IOUNT)
      ELSE
        TXML=XLOGL(IOUNT)
      ENDIF
      DO 1098 I=1,NIND
        MAX=1
        DO 1096 K=2,NG
1096      IF (W(I,K).GT.W(I,MAX)) MAX=K
        IDT(I)=MAX
        CKOUT=0.0
        DO 1097 K=1,NG
1097      CKOUT=CKOUT+W(I,K)
        IF (CKOUT.LT.0.0001) IDT(I)=0
1098  CONTINUE
      IF ((FLAGS(12).GT.0).AND.(FLAGS(9).EQ.1)) THEN
       CALL CAPART(NIND,NATT,NG,W,IDT,USA,XCC)     
        WRITE(FYLENO,*) 
        WRITE(FYLENO,*) '     *******************************'
       IF (FLAGS(8).EQ.4) THEN
        WRITE(FYLENO,*) '    FIT USING PARTIAL CLASSIFIED DATA'
       ELSE
        WRITE(FYLENO,*) '      FIT USING PARTIAL GROUPING '
       ENDIF
        WRITE(FYLENO,*) '     *******************************'
        WRITE(FYLENO,*) 
      WRITE(FYLENO,*)'  Implied grouping for all unclassified entities'
      WRITE(FYLENO,*)' (with component membership of classified'
      WRITE(FYLENO,*)'  entities as specified)'
        WRITE(FYLENO,1177) (IDT(III),III=1,NIND)
        WRITE (FYLENO,*)
1177    FORMAT (2X,10I4)
        CALL OUTESTIMATES(NIND,NATT,NG,NCOV,XMU,V,XVAR,DV,T,XUU)
        FLAGS(12)=0
        CALL ESTEP(NIND,NATT,NG,X,XMU,V,T,DEN,WL,W,XUU,USA,DV,
     &  XLOGL,IOUNT,XMAH,IER)
	IF (IER.GT.0) GOTO 9999 
        XTMP=1
        CALL MSTEP(NIND,NATT,NG,NCOV,X,XVAR,V,DV,
     &           XMU,WTOT,T,W,XUU,XMAH,XTMP,U,IER)
	IF (IER.GT.0) GOTO 9999 
        FLAGS(12)=1
        WRITE(FYLENO,*) 
        CALL OUTLOOP(NIND,NATT,NG,XMU,DV,T,NCOV,IOUNT,XLOGL,
     &	             W,IDT,X,DEN,USA,U)
      ELSEIF (FLAGS(9).EQ.1) THEN
         WRITE (FYLENO,1041) XLOGL(IOUNT)
1041     FORMAT(2X,'Final Log-Likelihood is ',F15.3)
         IF (FLAGS(11).EQ.0) THEN
           IF ((NG.NE.1).OR.(FLAGS(7).NE.0)) THEN
            CALL OUTLOOP(NIND,NATT,NG,XMU,DV,T,NCOV,IOUNT,XLOGL,
     &	                 W,IDT,X,DEN,USA,U)
           ENDIF
         ELSEIF (FLAGS(20).EQ.2) THEN
          DO 1172 III=1,NIND
           WRITE (29,1171) III,DEN(III,NG+1)
1171       FORMAT (2X,I6,2X,15F8.3)
1172      CONTINUE
         ENDIF
      ENDIF
9999   CONTINUE
      RETURN
      END

      SUBROUTINE CAPART(NIND,NATT,NG,W,IDT,USA,XCC)     
C     This subroutine determines  the partition of entities,
C     from the posterior probabilities W, into NG groups
      implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO,USA(MNIND)
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION W(MNIND,MAXNG),IDT(MNIND),
     &          XCC(MAXNG)
      DO 1110 K=1,NG
1110    XCC(K)=0.0
      DO 1140 I=1,NIND
        MAX=1
        DO 1120 K=2,NG
1120      IF (W(I,K).GT.W(I,MAX)) MAX=K
        XCC(MAX)=XCC(MAX)+W(I,MAX)
        IDT(I)=MAX
        CKOUT=0.0
        DO 1130 K=1,NG
1130      CKOUT=CKOUT+W(I,K)
        IF (CKOUT.LT.0.0001) IDT(I)=0
1140  CONTINUE
      RETURN
      END

      SUBROUTINE OUTLOOP(NIND,NATT,NG,XMU,DV,T,NCOV,IOUNT,XLOGL,
     &                   W,IDT,X,DEN,USA,U)     
c      This subroutine displays all the relevant information from the
c      EM algorithm applied to the best partition. 
      implicit double precision (a-h,o-z)
       INCLUDE 'EMMIX.max'
	      INTEGER FLAGS(40),FYLENO,USA(MNIND)
	      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION DV(MAXNG),T(MAXNG),W(MNIND,MAXNG),
     &          XCC(MAXNG),XLOGL(MITER),XMU(MAXNG,MNATT),
     &          IDT(MNIND),X(MNIND,MNATT),DEN(MNIND,MAXNG+1),
     &          U(MNIND,MNATT)
      IF (FLAGS(8).NE.4) THEN 
	WRITE (FYLENO,1105) IOUNT,XLOGL(IOUNT)
1105    FORMAT (/2X,'In loop ',I3,' log likelihood is ',F15.3)
      ENDIF
      CALL CAPART(NIND,NATT,NG,W,IDT,USA,XCC)     
      IF (FLAGS(6).EQ.1) THEN
        WRITE (*,*)
        WRITE (FYLENO,*) '  Density values of user defined model'
      WRITE(FYLENO,*)'Observation| component densities |',
     &' mixture log density '
        DO 1143 I=1,NIND
         IF (FLAGS(20).EQ.2) THEN
           WRITE (29,1141) I,(DEN(I,J),J=1,NG+1)
         ENDIF
1141       FORMAT (2X,I6,10G15.10)
           WRITE (FYLENO,1170) I,(DEN(I,J),J=1,NG+1)
1142       FORMAT (2X,I6,2X,10G15.5)
1143    CONTINUE
      ENDIF
        IF (FLAGS(12).EQ.1) THEN
      WRITE(FYLENO,*)'Observation | mixture log density',
     &' |Component1 Component2..etc..|',
     &'(specified component membership)'
        ELSEIF (FLAGS(7).GT.0) THEN
      WRITE (FYLENO,*) 'Observation | mixture log density|',
     &' Us Grp1 Grp2 ..etc...|'
        ELSE   
      WRITE (FYLENO,*) 'Observation | mixture log density|',
     &'Grp1 Grp2 ..etc...|'
        ENDIF
        DO 1160 I=1,NIND
          IF ((FLAGS(12).EQ.1).AND.(USA(I).GT.-1)) THEN
      WRITE (FYLENO,1170) I,DEN(I,NG+1),(W(I,K),K=1,NG),FLOAT(USA(I)) 
          ELSEIF (FLAGS(7).GT.0) THEN
            WRITE (FYLENO,1170) I,DEN(I,NG+1),(U(I,K),K=1,NG),
     &                 (W(I,K),K=1,NG)
          ELSE
            WRITE (FYLENO,1170) I,DEN(I,NG+1),(W(I,K),K=1,NG)
          ENDIF
1160    CONTINUE

c1170   FORMAT (2X,I6,2X,10G13.5,'*',I3)
1170   FORMAT (2X,I6,2X,G13.5,' ',10F7.4,'*',I3)
1171   FORMAT (2X,I6,2X,7G12.5)
      IF (FLAGS(12).EQ.0) THEN
      WRITE (FYLENO,1175) NG
1175  FORMAT (/2X,'Implied grouping of the entities into ',
     &     I3,' groups')
      ELSE
        WRITE (FYLENO,*) 
        WRITE (FYLENO,*)'  Implied grouping for all entities' 
        WRITE (FYLENO,*)'   (including classified entities)'
      ENDIF
      WRITE (FYLENO,1177) (IDT(I),I=1,NIND)
1177  FORMAT (2X,10I4)
      RETURN
      END



      SUBROUTINE LCAL(NIND,NATT,NG,X,W,XMU,WSUM,WUSUM,XVAR,U)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION X(MNIND,MNATT),W(MNIND,MAXNG),WUSUM(MAXNG),
     &          SUM(MAXNG,MNATT),XMU(MAXNG,MNATT),WSUM(MAXNG),
     &          XVAR(MAXNG,MNATT,MNATT),U(MNIND,MAXNG)
C     Compute estimate of mean, covariance matrix and mixing proportion
C     for each group
C     Compute new estimates of group means (XMU)
      DO 1310 K=1,NG
        DO 1310 J=1,NATT
          SUM(K,J)=0.0
          DO 1300 JJ=1,NIND
           SUM(K,J)=SUM(K,J)+X(JJ,J)*W(JJ,K)*U(JJ,K)
1300      CONTINUE
          XMU(K,J)=SUM(K,J)/WUSUM(K)
1310  CONTINUE
C     Compute new estimate of covariance matrix for each group
      DO 1360 K=1,NG
        DO 1320 J=1,NATT
          DO 1320 I=1,J
            XVAR(K,I,J)=0.0
1320    CONTINUE
        DO 1340 JJ=1,NIND
          DO 1330 J=1,NATT
            DO 1330 I=1,J
              XVAR(K,I,J)=XVAR(K,I,J)+(X(JJ,I)-XMU(K,I))
     &                  *(X(JJ,J)-XMU(K,J))*W(JJ,K)*U(JJ,K)
1330      CONTINUE
1340    CONTINUE
        DO 1350 J=1,NATT
          DO 1350 I=1,J
            XVAR(K,I,J)=XVAR(K,I,J)/WSUM(K)
            XVAR(K,J,I)=XVAR(K,I,J)
1350    CONTINUE
1360  CONTINUE
      RETURN
      END


      SUBROUTINE CRITERIA(NG,XLGL,NIND,NATT,NCOV,AIC,BIC,AWE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
C
C    This subroutine implements various criterion for determining the number of 
C    components. See Celeux and Soromenho "An Entropy Criterion for assessing
C    the Number of Clusters in a Mixture Model".                 
C                      Implemented by David Peel Sept 1994                     
C
       IF (NCOV.EQ.1) THEN
         VK =(NG-1) + NATT*NG + NATT*(NATT+1)/2
       ELSEIF (NCOV.EQ.2) THEN
         VK =(NG-1) + NATT*NG + NG*NATT*(NATT+1)/2
       ELSEIF (NCOV.EQ.3) THEN
         VK =(NG-1) + NATT*NG + NATT
       ELSEIF (NCOV.EQ.4) THEN
         VK =(NG-1) + NATT*NG + NG*NATT
       ENDIF
C      Calculate the value of the Akaike Information
            AIC=(-2.0)*XLGL+2.0*VK
C      Calculate the value of the Bayesian Information
            BIC=(-2.0)*XLGL+VK*log(FLOAT(NIND))
C      Calculate the value of the Approximate Weight of Evidence
            AWE=(-2.0)*XLGL+2.0*VK*(3.0/2.0+LOG(FLOAT(NIND)))
       RETURN
       END



      SUBROUTINE SYMINV(A,N,C,W,NULLTY,IFAULT,RMAX,CMY,TOL)
C
C       Algorithm AS7, Applied Statistics, Vol.17, 1968.
C
C       ARGUMENTS:-
C       A()     = input, the symmetric matrix to be inverted, stored in
C                 lower triangular form
C       N       = input, order of the matrix
C       C()     = output, the inverse of a (a generalized inverse if c is
C                 singular), also stored in lower triangular.
C                 c and a may occupy the same locations.
C       W()     = workspace, dimension at least n.
C       NULLTY  = output, the rank deficiency of a.
C       IFAULT  = output, error indicator
C                       = 1 if n < 1
C                       = 2 if a is not +ve semi-definite
C                       = 0 otherwise
C       RMAX    = output, approximate bound on the accuracy of the diagonal
C                 elements of c.  e.g. if rmax = 1.e-04 then the diagonal
C                 elements of c will be accurate to about 4 dec. digits.
C
C       Latest revision - 18 April 1981
 
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
        DIMENSION A(MNATT*MNATT),C(MNATT*MNATT),W(N),
     &          CMY(MNATT*(MNATT+1)/2)

        NROW=N
        IFAULT=1
        IF(NROW.LE.0) GO TO 100
        IFAULT=0
C       Cholesky factorization of A, result in C
        CALL CHOLA(A,NROW,C,NULLTY,IFAULT,RMAX,W,TOL)
        IF(IFAULT.NE.0) GO TO 100
C       Invert C & form the product (CINV)'*CINV, where CINV is the inverse
C       of C, row by row starting with the last row.
C       IROW = the row number, NDIAG = location of last element in the row.
        NN=NROW*(NROW+1)/2
        DO 200 IMY=1,NN
200       CMY(IMY)=C(IMY)
        IROW=NROW
        NDIAG=NN
10      IF(C(NDIAG).EQ.0.0) GO TO 60
        L=NDIAG
        DO 20 I=IROW,NROW
          W(I)=C(L)
          L=L+I
20      CONTINUE
        ICOL=NROW
        JCOL=NN
        MDIAG=NN
30      L=JCOL
        X=0.0
        IF(ICOL.EQ.IROW) X=1.0/W(IROW)
        K=NROW
40      IF(K.EQ.IROW) GO TO 50
        X=X-W(K)*C(L)
        K=K-1
        L=L-1
        IF(L.GT.MDIAG) L=L-K+1
        GO TO 40
50      C(L)=X/W(IROW)
        IF(ICOL.EQ.IROW) GO TO 80
        MDIAG=MDIAG-ICOL
        ICOL=ICOL-1
        JCOL=JCOL-1
        GO TO 30
60      L=NDIAG
        DO 70 J=IROW,NROW
          C(L)=0.0
          L=L+J
70      CONTINUE
80      NDIAG=NDIAG-IROW
        IROW=IROW-1
        IF(IROW.NE.0) GO TO 10
100     RETURN
        END


        SUBROUTINE CHOLA(A,N,U,NULLTY,IFAULT,RMAX,R,TOL)
C
C       Algorithm AS6, Applied Statistics, Vol.17, 1968, with
C       modifications by A.J.Miller
C
C       ARGUMENTS:-
C       A()     = input, a +ve definite matrix stored in lower-triangular
C                 form.
C       N       = input, the order of A
C       U()     = output, a lower triangular matrix such that U*U' = A.
C                 A & U may occupy the same locations.
C       NULLTY  = output, the rank deficiency of A.
C       IFAULT  = output, error indicator
C                       = 1 if N < 1
C                       = 2 if A is not +ve semi-definite
C                       = 0 otherwise
C       RMAX    = output, an estimate of the relative accuracy of the
C                 diagonal elements of U.
C       R()     = output, array containing bounds on the relative accuracy
C                 of each diagonal element of U.
C
C       Latest revision - 18 April 1981
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
        DIMENSION A(MNATT*MNATT),U(MNATT*MNATT),R(N)
 
C       ETA should be set equal to the smallest +ve value such that
C       1.0 + eta is calculated as being greater than 1.0 in the accuracy
C       being used.
C        DATA ETA/1.e-07/
        ETA=TOL
        IFAULT=1
        IF(N.LE.0) GO TO 100
        IFAULT=2
        NULLTY=0
        RMAX=ETA
        R(1)=ETA
        J=1
        K=0
C       Factorize column by column, ICOL = Column No.
        DO 80 ICOL=1,N
          L=0
C       IROW = row number within column ICOL
          DO 40 IROW=1,ICOL
            K=K+1
            W=A(K)
            IF(IROW.EQ.ICOL) RSQ=(W*ETA)**2
            M=J
            DO 10 I=1,IROW
              L=L+1
              IF(I.EQ.IROW) GO TO 20
              W=W-U(L)*U(M)
              IF(IROW.EQ.ICOL) RSQ=RSQ+(U(L)**2*R(I))**2
              M=M+1
10          CONTINUE
20          IF(IROW.EQ.ICOL) GO TO 50
            IF(U(L).EQ.0.0) GO TO 30
            U(K)=W/U(L)
            GO TO 40
30          U(K)=0.0
            IF(ABS(W).GT.ABS(RMAX*A(K))) GO TO 100
40        CONTINUE
C         End of row, estimate relative accuracy of diagonal element.
50        RSQ=SQRT(RSQ)
          IF(ABS(W).LE.5.*RSQ) GO TO 60
          IF(W.LT.0.0) GO TO 100
          U(K)=SQRT(W)
          R(I)=RSQ/W
          IF(R(I).GT.RMAX) RMAX=R(I)
          GO TO 70
60        U(K)=0.0
          NULLTY=NULLTY+1
70        J=J+ICOL
80      CONTINUE
        IFAULT=0.0
100     RETURN
        END

      SUBROUTINE SEM(NIND,NATT,NG,W)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      DIMENSION W(MNIND,MAXNG),XL(MAXNG),WSEM(MNIND,MAXNG)
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM,R
      DO 100 I=1,MNIND
       DO 100 K=1,MAXNG
	 WSEM(I,K)=0.0
100   CONTINUE

      DO 120 I=1,NIND
        XLIM=0.0
	R=RANDNUM()
        DO 110 K=1,NG-1
          XLIM=XLIM+W(I,K)
          XL(K)=XLIM
          IF (R.LE.XL(K)) THEN
            WSEM(I,K)=1
            GOTO 120
          ENDIF
110     CONTINUE 
        XL(NG)=1.0
	WSEM(I,NG)=1
120   CONTINUE
	   DO 945 II=1,NIND
	     DO 945 KK=1,NG
	       W(II,KK)=WSEM(II,KK)
945         CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
      INTEGER N,LR
      INCLUDE 'EMMIX.max'
      DOUBLE PRECISION DELTA
      DOUBLE PRECISION R(LR),DIAG(MAXNG),QTB(MAXNG),X(MAXNG),
     &                 WA1(MAXNG),WA2(MAXNG)
C     **********
C
C     SUBROUTINE DOGLEG
C
C     GIVEN AN M BY N MATRIX A, AN N BY N NONSINGULAR DIAGONAL
C     MATRIX D, AN M-VECTOR B, AND A POSITIVE NUMBER DELTA, THE
C     PROBLEM IS TO DETERMINE THE CONVEX COMBINATION X OF THE
C     GAUSS-NEWTON AND SCALED GRADIENT DIRECTIONS THAT MINIMIZES
C     (A*X - B) IN THE LEAST SQUARES SENSE, SUBJECT TO THE
C     RESTRICTION THAT THE EUCLIDEAN NORM OF D*X BE AT MOST DELTA.
C
C     THIS SUBROUTINE COMPLETES THE SOLUTION OF THE PROBLEM
C     IF IT IS PROVIDED WITH THE NECESSARY INFORMATION FROM THE
C     QR FACTORIZATION OF A. THAT IS, IF A = Q*R, WHERE Q HAS
C     ORTHOGONAL COLUMNS AND R IS AN UPPER TRIANGULAR MATRIX,
C     THEN DOGLEG EXPECTS THE FULL UPPER TRIANGLE OF R AND
C     THE FIRST N COMPONENTS OF (Q TRANSPOSE)*B.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE DOGLEG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE ORDER OF R.
C
C       R IS AN INPUT ARRAY OF LENGTH LR WHICH MUST CONTAIN THE UPPER
C         TRIANGULAR MATRIX R STORED BY ROWS.
C
C       LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
C         (N*(N+1))/2.
C
C       DIAG IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
C         DIAGONAL ELEMENTS OF THE MATRIX D.
C
C       QTB IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE FIRST
C         N ELEMENTS OF THE VECTOR (Q TRANSPOSE)*B.
C
C       DELTA IS A POSITIVE INPUT VARIABLE WHICH SPECIFIES AN UPPER
C         BOUND ON THE EUCLIDEAN NORM OF D*X.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE DESIRED
C         CONVEX COMBINATION OF THE GAUSS-NEWTON DIRECTION AND THE
C         SCALED GRADIENT DIRECTION.
C
C       WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N.
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... DPMPAR,ENORM
C
C       FORTRAN-SUPPLIED ... DABS,DMAX1,DMIN1,DSQRT
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,J,JJ,JP1,K,L
      DOUBLE PRECISION ALPHA,BNORM,EPSMCH,GNORM,ONE,QNORM,SGNORM,SUM,
     *                 TEMP,ZERO
      DOUBLE PRECISION DPMPAR,ENORM
      DATA ONE,ZERO /1.0D0,0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
      EPSMCH = DPMPAR(1)
C
C     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
C
      JJ = (N*(N + 1))/2 + 1
      DO 50 K = 1, N
         J = N - K + 1
         JP1 = J + 1
         JJ = JJ - K
         L = JJ + 1
         SUM = ZERO
         IF (N .LT. JP1) GO TO 20
         DO 10 I = JP1, N
            SUM = SUM + R(L)*X(I)
            L = L + 1
   10       CONTINUE
   20    CONTINUE
         TEMP = R(JJ)
         IF (TEMP .NE. ZERO) GO TO 40
         L = J
         DO 30 I = 1, J
            TEMP = DMAX1(TEMP,DABS(R(L)))
            L = L + N - I
   30       CONTINUE
         TEMP = EPSMCH*TEMP
         IF (TEMP .EQ. ZERO) TEMP = EPSMCH
   40    CONTINUE
         X(J) = (QTB(J) - SUM)/TEMP
   50    CONTINUE
C
C     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
C
      DO 60 J = 1, N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*X(J)
   60    CONTINUE
      QNORM = ENORM(N,WA2)
      IF (QNORM .LE. DELTA) GO TO 140
C
C     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
C     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
C
      L = 1
      DO 80 J = 1, N
         TEMP = QTB(J)
         DO 70 I = J, N
            WA1(I) = WA1(I) + R(L)*TEMP
            L = L + 1
   70       CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   80    CONTINUE
C
C     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
C     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
C
      GNORM = ENORM(N,WA1)
      SGNORM = ZERO
      ALPHA = DELTA/QNORM
      IF (GNORM .EQ. ZERO) GO TO 120
C
C     CALCULATE THE POINT ALONG THE SCALED GRADIENT
C     AT WHICH THE QUADRATIC IS MINIMIZED.
C
      DO 90 J = 1, N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   90    CONTINUE
      L = 1
      DO 110 J = 1, N
         SUM = ZERO
         DO 100 I = J, N
            SUM = SUM + R(L)*WA1(I)
            L = L + 1
  100       CONTINUE
         WA2(J) = SUM
  110    CONTINUE
      TEMP = ENORM(N,WA2)
      SGNORM = (GNORM/TEMP)/TEMP
C
C     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
C
      ALPHA = ZERO
      IF (SGNORM .GE. DELTA) GO TO 120
C
C     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
C     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
C     AT WHICH THE QUADRATIC IS MINIMIZED.
C
      BNORM = ENORM(N,QTB)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2
     *       + DSQRT((TEMP-(DELTA/QNORM))**2
     *               +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP
  120 CONTINUE
C
C     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
C     DIRECTION AND THE SCALED GRADIENT DIRECTION.
C
      TEMP = (ONE - ALPHA)*DMIN1(SGNORM,DELTA)
      DO 130 J = 1, N
         X(J) = TEMP*WA1(J) + ALPHA*X(J)
  130    CONTINUE
  140 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE DOGLEG.
C
      END
      DOUBLE PRECISION FUNCTION DPMPAR(I)
      INTEGER I
C     **********
C
C     FUNCTION DPMPAR
C
C     THIS FUNCTION PROVIDES DOUBLE PRECISION MACHINE PARAMETERS
C     WHEN THE APPROPRIATE SET OF DATA STATEMENTS IS ACTIVATED (BY
C     REMOVING THE C FROM COLUMN 1) AND ALL OTHER DATA STATEMENTS ARE
C     RENDERED INACTIVE. MOST OF THE PARAMETER VALUES WERE OBTAINED
C     FROM THE CORRESPONDING BELL LABORATORIES PORT LIBRARY FUNCTION.
C
C     THE FUNCTION STATEMENT IS
C
C       DOUBLE PRECISION FUNCTION DPMPAR(I)
C
C     WHERE
C
C       I IS AN INTEGER INPUT VARIABLE SET TO 1, 2, OR 3 WHICH
C         SELECTS THE DESIRED MACHINE PARAMETER. IF THE MACHINE HAS
C         T BASE B DIGITS AND ITS SMALLEST AND LARGEST EXPONENTS ARE
C         EMIN AND EMAX, RESPECTIVELY, THEN THESE PARAMETERS ARE
C
C         DPMPAR(1) = B**(1 - T), THE MACHINE PRECISION,
C
C         DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
C
C         DPMPAR(3) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER MCHEPS(4)
      INTEGER MINMAG(4)
      INTEGER MAXMAG(4)
      DOUBLE PRECISION DMACH(3)
      EQUIVALENCE (DMACH(1),MCHEPS(1))
      EQUIVALENCE (DMACH(2),MINMAG(1))
      EQUIVALENCE (DMACH(3),MAXMAG(1))
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE AMDAHL 470/V6, THE ICL 2900, THE ITEL AS/6,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C     DATA MCHEPS(1),MCHEPS(2) / Z34100000, Z00000000 /
C     DATA MINMAG(1),MINMAG(2) / Z00100000, Z00000000 /
C     DATA MAXMAG(1),MAXMAG(2) / Z7FFFFFFF, ZFFFFFFFF /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
C
C     DATA MCHEPS(1),MCHEPS(2) / O606400000000, O000000000000 /
C     DATA MINMAG(1),MINMAG(2) / O402400000000, O000000000000 /
C     DATA MAXMAG(1),MAXMAG(2) / O376777777777, O777777777777 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA MCHEPS(1) / 15614000000000000000B /
C     DATA MCHEPS(2) / 15010000000000000000B /
C
C     DATA MINMAG(1) / 00604000000000000000B /
C     DATA MINMAG(2) / 00000000000000000000B /
C
C     DATA MAXMAG(1) / 37767777777777777777B /
C     DATA MAXMAG(2) / 37167777777777777777B /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA MCHEPS(1),MCHEPS(2) / "114400000000, "000000000000 /
C     DATA MINMAG(1),MINMAG(2) / "033400000000, "000000000000 /
C     DATA MAXMAG(1),MAXMAG(2) / "377777777777, "344777777777 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA MCHEPS(1),MCHEPS(2) / "104400000000, "000000000000 /
C     DATA MINMAG(1),MINMAG(2) / "000400000000, "000000000000 /
C     DATA MAXMAG(1),MAXMAG(2) / "377777777777, "377777777777 /
C
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA MCHEPS(1),MCHEPS(2) /  620756992,           0 /
C     DATA MINMAG(1),MINMAG(2) /    8388608,           0 /
C     DATA MAXMAG(1),MAXMAG(2) / 2147483647,          -1 /
C
C     DATA MCHEPS(1),MCHEPS(2) / O04500000000, O00000000000 /
C     DATA MINMAG(1),MINMAG(2) / O00040000000, O00000000000 /
C     DATA MAXMAG(1),MAXMAG(2) / O17777777777, O37777777777 /
C
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA MCHEPS(1),MCHEPS(2) /   9472,      0 /
C     DATA MCHEPS(3),MCHEPS(4) /      0,      0 /
C
C     DATA MINMAG(1),MINMAG(2) /    128,      0 /
C     DATA MINMAG(3),MINMAG(4) /      0,      0 /
C
C     DATA MAXMAG(1),MAXMAG(2) /  32767,     -1 /
C     DATA MAXMAG(3),MAXMAG(4) /     -1,     -1 /
C
C     DATA MCHEPS(1),MCHEPS(2) / O022400, O000000 /
C     DATA MCHEPS(3),MCHEPS(4) / O000000, O000000 /
C
C     DATA MINMAG(1),MINMAG(2) / O000200, O000000 /
C     DATA MINMAG(3),MINMAG(4) / O000000, O000000 /
C
C     DATA MAXMAG(1),MAXMAG(2) / O077777, O177777 /
C     DATA MAXMAG(3),MAXMAG(4) / O177777, O177777 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C     DATA MCHEPS(1) / O1451000000000000 /
C     DATA MCHEPS(2) / O0000000000000000 /
C
C     DATA MINMAG(1) / O1771000000000000 /
C     DATA MINMAG(2) / O7770000000000000 /
C
C     DATA MAXMAG(1) / O0777777777777777 /
C     DATA MAXMAG(2) / O7777777777777777 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C     DATA MCHEPS(1) / O1451000000000000 /
C     DATA MCHEPS(2) / O0000000000000000 /
C
C     DATA MINMAG(1) / O1771000000000000 /
C     DATA MINMAG(2) / O0000000000000000 /
C
C     DATA MAXMAG(1) / O0777777777777777 /
C     DATA MAXMAG(2) / O0007777777777777 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA MCHEPS(1) / ZCC6800000 /
C     DATA MCHEPS(2) / Z000000000 /
C
C     DATA MINMAG(1) / ZC00800000 /
C     DATA MINMAG(2) / Z000000000 /
C
C     DATA MAXMAG(1) / ZDFFFFFFFF /
C     DATA MAXMAG(2) / ZFFFFFFFFF /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     DATA MCHEPS(1),MCHEPS(2) / O170640000000, O000000000000 /
C     DATA MINMAG(1),MINMAG(2) / O000040000000, O000000000000 /
C     DATA MAXMAG(1),MAXMAG(2) / O377777777777, O777777777777 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC DMACH(3)
C
C     DATA MINMAG/20K,3*0/,MAXMAG/77777K,3*177777K/
C     DATA MCHEPS/32020K,3*0/
C
C     MACHINE CONSTANTS FOR THE HARRIS 220.
C
C     DATA MCHEPS(1),MCHEPS(2) / '20000000, '00000334 /
C     DATA MINMAG(1),MINMAG(2) / '20000000, '00000201 /
C     DATA MAXMAG(1),MAXMAG(2) / '37777777, '37777577 /
C
C     MACHINE CONSTANTS FOR THE CRAY-1.
C
C     DATA MCHEPS(1) / 0376424000000000000000B /
C     DATA MCHEPS(2) / 0000000000000000000000B /
C
C     DATA MINMAG(1) / 0200034000000000000000B /
C     DATA MINMAG(2) / 0000000000000000000000B /
C
C     DATA MAXMAG(1) / 0577777777777777777777B /
C     DATA MAXMAG(2) / 0000007777777777777776B /
C
C     MACHINE CONSTANTS FOR THE PRIME 400.
C
C     DATA MCHEPS(1),MCHEPS(2) / :10000000000, :00000000123 /
C     DATA MINMAG(1),MINMAG(2) / :10000000000, :00000100000 /
C     DATA MAXMAG(1),MAXMAG(2) / :17777777777, :37777677776 /
C
C     MACHINE CONSTANTS FOR THE VAX-11.
C
C     DATA MCHEPS(1),MCHEPS(2) /   9472,  0 /
C     DATA MINMAG(1),MINMAG(2) /    128,  0 /
C     DATA MAXMAG(1),MAXMAG(2) / -32769, -1 /
C
C     DPMPAR = DMACH(I)
      if(i.eq.1) then
	DPMPAR=    2.2204460492503D-16
      else if (i.eq.2) then
	DPMPAR=    2.2250738585072D-308
      else
	DPMPAR=    8.9884656743116D+307
      endif
      RETURN
C
C     LAST CARD OF FUNCTION DPMPAR.
C
      END
      DOUBLE PRECISION FUNCTION ENORM(N,X)
      INCLUDE 'EMMIX.max'
      INTEGER N
      DOUBLE PRECISION X(MAXNG)
C     **********
C
C     FUNCTION ENORM
C
C     GIVEN AN N-VECTOR X, THIS FUNCTION CALCULATES THE
C     EUCLIDEAN NORM OF X.
C
C     THE EUCLIDEAN NORM IS COMPUTED BY ACCUMULATING THE SUM OF
C     SQUARES IN THREE DIFFERENT SUMS. THE SUMS OF SQUARES FOR THE
C     SMALL AND LARGE COMPONENTS ARE SCALED SO THAT NO OVERFLOWS
C     OCCUR. NON-DESTRUCTIVE UNDERFLOWS ARE PERMITTED. UNDERFLOWS
C     AND OVERFLOWS DO NOT OCCUR IN THE COMPUTATION OF THE UNSCALED
C     SUM OF SQUARES FOR THE INTERMEDIATE COMPONENTS.
C     THE DEFINITIONS OF SMALL, INTERMEDIATE AND LARGE COMPONENTS
C     DEPEND ON TWO CONSTANTS, RDWARF AND RGIANT. THE MAIN
C     RESTRICTIONS ON THESE CONSTANTS ARE THAT RDWARF**2 NOT
C     UNDERFLOW AND RGIANT**2 NOT OVERFLOW. THE CONSTANTS
C     GIVEN HERE ARE SUITABLE FOR EVERY KNOWN COMPUTER.
C
C     THE FUNCTION STATEMENT IS
C
C       DOUBLE PRECISION FUNCTION ENORM(N,X)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DABS,DSQRT
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I
      DOUBLE PRECISION AGIANT,FLOATN,ONE,RDWARF,RGIANT,S1,S2,S3,XABS,
     *                 X1MAX,X3MAX,ZERO
      DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = DABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
C
C              SUM FOR LARGE COMPONENTS.
C
               IF (XABS .LE. X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
C
C              SUM FOR SMALL COMPONENTS.
C
               IF (XABS .LE. X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
C
C           SUM FOR INTERMEDIATE COMPONENTS.
C
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
C
C     CALCULATION OF NORM.
C
      IF (S1 .EQ. ZERO) GO TO 100
         ENORM = X1MAX*DSQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (S2 .EQ. ZERO) GO TO 110
            IF (S2 .GE. X3MAX)
     *         ENORM = DSQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 .LT. X3MAX)
     *         ENORM = DSQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            ENORM = X3MAX*DSQRT(S3)
  120    CONTINUE
  130 CONTINUE
      RETURN
C
C     LAST CARD OF FUNCTION ENORM.
C
      END
      SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
     *                  WA1,WA2)
      INCLUDE 'EMMIX.max'
      INTEGER N,LDFJAC,IFLAG,ML,MU
      DOUBLE PRECISION EPSFCN
      DOUBLE PRECISION X(MAXNG),FVEC(*),FJAC(LDFJAC,MAXNG),WA1(MAXNG),
     &                 WA2(MAXNG)
C     **********
C
C     SUBROUTINE FDJAC1
C
C     THIS SUBROUTINE COMPUTES A FORWARD-DIFFERENCE APPROXIMATION
C     TO THE N BY N JACOBIAN MATRIX ASSOCIATED WITH A SPECIFIED
C     PROBLEM OF N FUNCTIONS IN N VARIABLES. IF THE JACOBIAN HAS
C     A BANDED FORM, THEN FUNCTION EVALUATIONS ARE SAVED BY ONLY
C     APPROXIMATING THE NONZERO TERMS.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
C                         WA1,WA2)
C
C     WHERE
C
C       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
C         CALCULATES THE FUNCTIONS. FCN MUST BE DECLARED
C         IN AN EXTERNAL STATEMENT IN THE USER CALLING
C         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(MAXNG),FVEC(*)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ----------
C         RETURN
C         END
C
C         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
C         THE USER WANTS TO TERMINATE EXECUTION OF FDJAC1.
C         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF FUNCTIONS AND VARIABLES.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       FVEC IS AN INPUT ARRAY OF LENGTH N WHICH MUST CONTAIN THE
C         FUNCTIONS EVALUATED AT X.
C
C       FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE
C         APPROXIMATION TO THE JACOBIAN MATRIX EVALUATED AT X.
C
C       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
C
C       IFLAG IS AN INTEGER VARIABLE WHICH CAN BE USED TO TERMINATE
C         THE EXECUTION OF FDJAC1. SEE DESCRIPTION OF FCN.
C
C       ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
C         THE NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE
C         JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
C         ML TO AT LEAST N - 1.
C
C       EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
C         STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
C         APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
C         FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
C         THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
C         ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
C         PRECISION.
C
C       MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
C         THE NUMBER OF SUPERDIAGONALS WITHIN THE BAND OF THE
C         JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
C         MU TO AT LEAST N - 1.
C
C       WA1 AND WA2 ARE WORK ARRAYS OF LENGTH N. IF ML + MU + 1 IS AT
C         LEAST N, THEN THE JACOBIAN IS CONSIDERED DENSE, AND WA2 IS
C         NOT REFERENCED.
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... DPMPAR
C
C       FORTRAN-SUPPLIED ... DABS,DMAX1,DSQRT
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,J,K,MSUM
      DOUBLE PRECISION EPS,EPSMCH,H,TEMP,ZERO
      DOUBLE PRECISION DPMPAR
      DATA ZERO /0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
      EPSMCH = DPMPAR(1)
C
      EPS = DSQRT(DMAX1(EPSFCN,EPSMCH))
      MSUM = ML + MU + 1
      IF (MSUM .LT. N) GO TO 40
C
C        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
C
         DO 20 J = 1, N
            TEMP = X(J)
            H = EPS*DABS(TEMP)
            IF (H .EQ. ZERO) H = EPS
            X(J) = TEMP + H
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG .LT. 0) GO TO 30
            X(J) = TEMP
            DO 10 I = 1, N
               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
         GO TO 110
   40 CONTINUE
C
C        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
C
         DO 90 K = 1, MSUM
            DO 60 J = K, N, MSUM
               WA2(J) = X(J)
               H = EPS*DABS(WA2(J))
               IF (H .EQ. ZERO) H = EPS
               X(J) = WA2(J) + H
   60          CONTINUE
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG .LT. 0) GO TO 100
            DO 80 J = K, N, MSUM
               X(J) = WA2(J)
               H = EPS*DABS(WA2(J))
               IF (H .EQ. ZERO) H = EPS
               DO 70 I = 1, N
                  FJAC(I,J) = ZERO
                  IF (I .GE. J - MU .AND. I .LE. J + ML)
     *               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE FDJAC1.
C
      END
      SUBROUTINE HYBRD(FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,DIAG,
     *                 MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,LDFJAC,R,LR,
     *                 QTF,WA1,WA2,WA3,WA4)
      INCLUDE 'EMMIX.max'
      INTEGER N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR
      DOUBLE PRECISION XTOL,EPSFCN,FACTOR
      DOUBLE PRECISION X(MAXNG),FVEC(*),DIAG(MAXNG),FJAC(LDFJAC,MAXNG),
     &               R(LR),QTF(MAXNG),WA1(MAXNG),WA2(MAXNG),WA3(MAXNG),
     &                 WA4(MAXNG)
      EXTERNAL FCN
C     **********
C
C     SUBROUTINE HYBRD
C
C     THE PURPOSE OF HYBRD IS TO FIND A ZERO OF A SYSTEM OF
C     N NONLINEAR FUNCTIONS IN N VARIABLES BY A MODIFICATION
C     OF THE POWELL HYBRID METHOD. THE USER MUST PROVIDE A
C     SUBROUTINE WHICH CALCULATES THE FUNCTIONS. THE JACOBIAN IS
C     THEN CALCULATED BY A FORWARD-DIFFERENCE APPROXIMATION.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE HYBRD(FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,
C                        DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,FJAC,
C                        LDFJAC,R,LR,QTF,WA1,WA2,WA3,WA4)
C
C     WHERE
C
C       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
C         CALCULATES THE FUNCTIONS. FCN MUST BE DECLARED
C         IN AN EXTERNAL STATEMENT IN THE USER CALLING
C         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(MAXNG),FVEC(*)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ---------
C         RETURN
C         END
C
C         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
C         THE USER WANTS TO TERMINATE EXECUTION OF HYBRD.
C         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF FUNCTIONS AND VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN
C         AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X
C         CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR.
C
C       FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
C         THE FUNCTIONS EVALUATED AT THE OUTPUT X.
C
C       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION
C         OCCURS WHEN THE RELATIVE ERROR BETWEEN TWO CONSECUTIVE
C         ITERATES IS AT MOST XTOL.
C
C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
C         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV
C         BY THE END OF AN ITERATION.
C
C       ML IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
C         THE NUMBER OF SUBDIAGONALS WITHIN THE BAND OF THE
C         JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
C         ML TO AT LEAST N - 1.
C
C       MU IS A NONNEGATIVE INTEGER INPUT VARIABLE WHICH SPECIFIES
C         THE NUMBER OF SUPERDIAGONALS WITHIN THE BAND OF THE
C         JACOBIAN MATRIX. IF THE JACOBIAN IS NOT BANDED, SET
C         MU TO AT LEAST N - 1.
C
C       EPSFCN IS AN INPUT VARIABLE USED IN DETERMINING A SUITABLE
C         STEP LENGTH FOR THE FORWARD-DIFFERENCE APPROXIMATION. THIS
C         APPROXIMATION ASSUMES THAT THE RELATIVE ERRORS IN THE
C         FUNCTIONS ARE OF THE ORDER OF EPSFCN. IF EPSFCN IS LESS
C         THAN THE MACHINE PRECISION, IT IS ASSUMED THAT THE RELATIVE
C         ERRORS IN THE FUNCTIONS ARE OF THE ORDER OF THE MACHINE
C         PRECISION.
C
C       DIAG IS AN ARRAY OF LENGTH N. IF MODE = 1 (SEE
C         BELOW), DIAG IS INTERNALLY SET. IF MODE = 2, DIAG
C         MUST CONTAIN POSITIVE ENTRIES THAT SERVE AS
C         MULTIPLICATIVE SCALE FACTORS FOR THE VARIABLES.
C
C       MODE IS AN INTEGER INPUT VARIABLE. IF MODE = 1, THE
C         VARIABLES WILL BE SCALED INTERNALLY. IF MODE = 2,
C         THE SCALING IS SPECIFIED BY THE INPUT DIAG. OTHER
C         VALUES OF MODE ARE EQUIVALENT TO MODE = 1.
C
C       FACTOR IS A POSITIVE INPUT VARIABLE USED IN DETERMINING THE
C         INITIAL STEP BOUND. THIS BOUND IS SET TO THE PRODUCT OF
C         FACTOR AND THE EUCLIDEAN NORM OF DIAG*X IF NONZERO, OR ELSE
C         TO FACTOR ITSELF. IN MOST CASES FACTOR SHOULD LIE IN THE
C         INTERVAL (.1,100.). 100. IS A GENERALLY RECOMMENDED VALUE.
C
C       NPRINT IS AN INTEGER INPUT VARIABLE THAT ENABLES CONTROLLED
C         PRINTING OF ITERATES IF IT IS POSITIVE. IN THIS CASE,
C         FCN IS CALLED WITH IFLAG = 0 AT THE BEGINNING OF THE FIRST
C         ITERATION AND EVERY NPRINT ITERATIONS THEREAFTER AND
C         IMMEDIATELY PRIOR TO RETURN, WITH X AND FVEC AVAILABLE
C         FOR PRINTING. IF NPRINT IS NOT POSITIVE, NO SPECIAL CALLS
C         OF FCN WITH IFLAG = 0 ARE MADE.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
C         TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
C         VALUE OF IFLAG. SEE DESCRIPTION OF FCN. OTHERWISE,
C         INFO IS SET AS FOLLOWS.
C
C         INFO = 0   IMPROPER INPUT PARAMETERS.
C
C         INFO = 1   RELATIVE ERROR BETWEEN TWO CONSECUTIVE ITERATES
C                    IS AT MOST XTOL.
C
C         INFO = 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED
C                    MAXFEV.
C
C         INFO = 3   XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
C                    THE APPROXIMATE SOLUTION X IS POSSIBLE.
C
C         INFO = 4   ITERATION IS NOT MAKING GOOD PROGRESS, AS
C                    MEASURED BY THE IMPROVEMENT FROM THE LAST
C                    FIVE JACOBIAN EVALUATIONS.
C
C         INFO = 5   ITERATION IS NOT MAKING GOOD PROGRESS, AS
C                    MEASURED BY THE IMPROVEMENT FROM THE LAST
C                    TEN ITERATIONS.
C
C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
C         CALLS TO FCN.
C
C       FJAC IS AN OUTPUT N BY N ARRAY WHICH CONTAINS THE
C         ORTHOGONAL MATRIX Q PRODUCED BY THE QR FACTORIZATION
C         OF THE FINAL APPROXIMATE JACOBIAN.
C
C       LDFJAC IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN N
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY FJAC.
C
C       R IS AN OUTPUT ARRAY OF LENGTH LR WHICH CONTAINS THE
C         UPPER TRIANGULAR MATRIX PRODUCED BY THE QR FACTORIZATION
C         OF THE FINAL APPROXIMATE JACOBIAN, STORED ROWWISE.
C
C       LR IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
C         (N*(N+1))/2.
C
C       QTF IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
C         THE VECTOR (Q TRANSPOSE)*FVEC.
C
C       WA1, WA2, WA3, AND WA4 ARE WORK ARRAYS OF LENGTH N.
C
C     SUBPROGRAMS CALLED
C
C       USER-SUPPLIED ...... FCN
C
C       MINPACK-SUPPLIED ... DOGLEG,DPMPAR,ENORM,FDJAC1,
C                            QFORM,QRFAC,R1MPYQ,R1UPDT
C
C       FORTRAN-SUPPLIED ... DABS,DMAX1,DMIN1,MIN0,MOD
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IFLAG,ITER,J,JM1,L,MSUM,NCFAIL,NCSUC,NSLOW1,NSLOW2
      INTEGER IWA(1)
      LOGICAL JEVAL,SING
      DOUBLE PRECISION ACTRED,DELTA,EPSMCH,FNORM,FNORM1,ONE,PNORM,
     *                 PRERED,P1,P5,P001,P0001,RATIO,SUM,TEMP,XNORM,
     *                 ZERO
      DOUBLE PRECISION DPMPAR,ENORM
      DATA ONE,P1,P5,P001,P0001,ZERO
     *     /1.0D0,1.0D-1,5.0D-1,1.0D-3,1.0D-4,0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
      EPSMCH = DPMPAR(1)
C
      INFO = 0
      IFLAG = 0
      NFEV = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. XTOL .LT. ZERO .OR. MAXFEV .LE. 0
     *    .OR. ML .LT. 0 .OR. MU .LT. 0 .OR. FACTOR .LE. ZERO
     *    .OR. LDFJAC .LT. N .OR. LR .LT. (N*(N + 1))/2) GO TO 300
      IF (MODE .NE. 2) GO TO 20
      DO 10 J = 1, N
         IF (DIAG(J) .LE. ZERO) GO TO 300
   10    CONTINUE
   20 CONTINUE
C
C     EVALUATE THE FUNCTION AT THE STARTING POINT
C     AND CALCULATE ITS NORM.
C
      IFLAG = 1
      CALL FCN(N,X,FVEC,IFLAG)
      NFEV = 1
      IF (IFLAG .LT. 0) GO TO 300
      FNORM = ENORM(N,FVEC)
C
C     DETERMINE THE NUMBER OF CALLS TO FCN NEEDED TO COMPUTE
C     THE JACOBIAN MATRIX.
C
      MSUM = MIN0(ML+MU+1,N)
C
C     INITIALIZE ITERATION COUNTER AND MONITORS.
C
      ITER = 1
      NCSUC = 0
      NCFAIL = 0
      NSLOW1 = 0
      NSLOW2 = 0
C
C     BEGINNING OF THE OUTER LOOP.
C
   30 CONTINUE
         JEVAL = .TRUE.
C
C        CALCULATE THE JACOBIAN MATRIX.
C
         IFLAG = 2
         CALL FDJAC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,WA1,
     *               WA2)
         NFEV = NFEV + MSUM
         IF (IFLAG .LT. 0) GO TO 300
C
C        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
         CALL QRFAC(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)
C
C        ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
C        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
C
         IF (ITER .NE. 1) GO TO 70
         IF (MODE .EQ. 2) GO TO 50
         DO 40 J = 1, N
            DIAG(J) = WA2(J)
            IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   40       CONTINUE
   50    CONTINUE
C
C        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
C        AND INITIALIZE THE STEP BOUND DELTA.
C
         DO 60 J = 1, N
            WA3(J) = DIAG(J)*X(J)
   60       CONTINUE
         XNORM = ENORM(N,WA3)
         DELTA = FACTOR*XNORM
         IF (DELTA .EQ. ZERO) DELTA = FACTOR
   70    CONTINUE
C
C        FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
C
         DO 80 I = 1, N
            QTF(I) = FVEC(I)
   80       CONTINUE
         DO 120 J = 1, N
            IF (FJAC(J,J) .EQ. ZERO) GO TO 110
            SUM = ZERO
            DO 90 I = J, N
               SUM = SUM + FJAC(I,J)*QTF(I)
   90          CONTINUE
            TEMP = (-SUM)/FJAC(J,J)
            DO 100 I = J, N
               QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  100          CONTINUE
  110       CONTINUE
  120       CONTINUE
C
C   COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
C
         SING = .FALSE.
         DO 150 J = 1, N
            L = J
            JM1 = J - 1
            IF (JM1 .LT. 1) GO TO 140
            DO 130 I = 1, JM1
               R(L) = FJAC(I,J)
               L = L + N - I
  130          CONTINUE
  140       CONTINUE
            R(L) = WA1(J)
            IF (WA1(J) .EQ. ZERO) SING = .TRUE.
  150       CONTINUE
C
C        ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
C
         CALL QFORM(N,N,FJAC,LDFJAC,WA1)
C
C        RESCALE IF NECESSARY.
C
         IF (MODE .EQ. 2) GO TO 170
         DO 160 J = 1, N
            DIAG(J) = DMAX1(DIAG(J),WA2(J))
  160       CONTINUE
  170    CONTINUE
C
C        BEGINNING OF THE INNER LOOP.
C
  180    CONTINUE
C
C           IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
C
            IF (NPRINT .LE. 0) GO TO 190
            IFLAG = 0
            IF (MOD(ITER-1,NPRINT) .EQ. 0) CALL FCN(N,X,FVEC,IFLAG)
            IF (IFLAG .LT. 0) GO TO 300
  190       CONTINUE
C
C           DETERMINE THE DIRECTION P.
C
            CALL DOGLEG(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
C
C           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
C
            DO 200 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
  200          CONTINUE
            PNORM = ENORM(N,WA3)
C
C           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
            IF (ITER .EQ. 1) DELTA = DMIN1(DELTA,PNORM)
C
C           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
C
            IFLAG = 1
            CALL FCN(N,WA2,WA4,IFLAG)
            NFEV = NFEV + 1
            IF (IFLAG .LT. 0) GO TO 300
            FNORM1 = ENORM(N,WA4)
C
C           COMPUTE THE SCALED ACTUAL REDUCTION.
C
            ACTRED = -ONE
            IF (FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C           COMPUTE THE SCALED PREDICTED REDUCTION.
C
            L = 1
            DO 220 I = 1, N
               SUM = ZERO
               DO 210 J = I, N
                  SUM = SUM + R(L)*WA1(J)
                  L = L + 1
  210             CONTINUE
               WA3(I) = QTF(I) + SUM
  220          CONTINUE
            TEMP = ENORM(N,WA3)
            PRERED = ZERO
            IF (TEMP .LT. FNORM) PRERED = ONE - (TEMP/FNORM)**2
C
C           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
C           REDUCTION.
C
            RATIO = ZERO
            IF (PRERED .GT. ZERO) RATIO = ACTRED/PRERED
C
C           UPDATE THE STEP BOUND.
C
            IF (RATIO .GE. P1) GO TO 230
               NCSUC = 0
               NCFAIL = NCFAIL + 1
               DELTA = P5*DELTA
               GO TO 240
  230       CONTINUE
               NCFAIL = 0
               NCSUC = NCSUC + 1
               IF (RATIO .GE. P5 .OR. NCSUC .GT. 1)
     *            DELTA = DMAX1(DELTA,PNORM/P5)
               IF (DABS(RATIO-ONE) .LE. P1) DELTA = PNORM/P5
  240       CONTINUE
C
C           TEST FOR SUCCESSFUL ITERATION.
C
            IF (RATIO .LT. P0001) GO TO 260
C
C           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
C
            DO 250 J = 1, N
               X(J) = WA2(J)
               WA2(J) = DIAG(J)*X(J)
               FVEC(J) = WA4(J)
  250          CONTINUE
            XNORM = ENORM(N,WA2)
            FNORM = FNORM1
            ITER = ITER + 1
  260       CONTINUE
C
C           DETERMINE THE PROGRESS OF THE ITERATION.
C
            NSLOW1 = NSLOW1 + 1
            IF (ACTRED .GE. P001) NSLOW1 = 0
            IF (JEVAL) NSLOW2 = NSLOW2 + 1
            IF (ACTRED .GE. P1) NSLOW2 = 0
C
C           TEST FOR CONVERGENCE.
C
            IF (DELTA .LE. XTOL*XNORM .OR. FNORM .EQ. ZERO) INFO = 1
            IF (INFO .NE. 0) GO TO 300
C
C           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
C
            IF (NFEV .GE. MAXFEV) INFO = 2
            IF (P1*DMAX1(P1*DELTA,PNORM) .LE. EPSMCH*XNORM) INFO = 3
            IF (NSLOW2 .EQ. 5) INFO = 4
            IF (NSLOW1 .EQ. 10) INFO = 5
            IF (INFO .NE. 0) GO TO 300
C
C           CRITERION FOR RECALCULATING JACOBIAN APPROXIMATION
C           BY FORWARD DIFFERENCES.
C
            IF (NCFAIL .EQ. 2) GO TO 290
C
C           CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
C           AND UPDATE QTF IF NECESSARY.
C
            DO 280 J = 1, N
               SUM = ZERO
               DO 270 I = 1, N
                  SUM = SUM + FJAC(I,J)*WA4(I)
  270             CONTINUE
               WA2(J) = (SUM - WA3(J))/PNORM
               WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
               IF (RATIO .GE. P0001) QTF(J) = SUM
  280          CONTINUE
C
C           COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
C
            CALL R1UPDT(N,N,R,LR,WA1,WA2,WA3,SING)
            CALL R1MPYQ(N,N,FJAC,LDFJAC,WA2,WA3)
            CALL R1MPYQ(1,N,QTF,1,WA2,WA3)
C
C           END OF THE INNER LOOP.
C
            JEVAL = .FALSE.
            GO TO 180
  290    CONTINUE
C
C        END OF THE OUTER LOOP.
C
         GO TO 30
  300 CONTINUE
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL FCN(N,X,FVEC,IFLAG)
      RETURN
C
C     LAST CARD OF SUBROUTINE HYBRD.
C
      END
      SUBROUTINE HYBRD1(FCN,N,X,FVEC,TOL,INFO,WA,LWA)
      INCLUDE 'EMMIX.max'
      INTEGER N,INFO,LWA
      DOUBLE PRECISION TOL
      DOUBLE PRECISION X(MAXNG),FVEC(*),WA((MAXNG*(3*MAXNG+13))/2)
      EXTERNAL FCN
C     **********
C
C     SUBROUTINE HYBRD1
C
C     THE PURPOSE OF HYBRD1 IS TO FIND A ZERO OF A SYSTEM OF
C     N NONLINEAR FUNCTIONS IN N VARIABLES BY A MODIFICATION
C     OF THE POWELL HYBRID METHOD. THIS IS DONE BY USING THE
C     MORE GENERAL NONLINEAR EQUATION SOLVER HYBRD. THE USER
C     MUST PROVIDE A SUBROUTINE WHICH CALCULATES THE FUNCTIONS.
C     THE JACOBIAN IS THEN CALCULATED BY A FORWARD-DIFFERENCE
C     APPROXIMATION.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE HYBRD1(FCN,N,X,FVEC,TOL,INFO,WA,LWA)
C
C     WHERE
C
C       FCN IS THE NAME OF THE USER-SUPPLIED SUBROUTINE WHICH
C         CALCULATES THE FUNCTIONS. FCN MUST BE DECLARED
C         IN AN EXTERNAL STATEMENT IN THE USER CALLING
C         PROGRAM, AND SHOULD BE WRITTEN AS FOLLOWS.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(MAXNG),FVEC(*)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ---------
C         RETURN
C         END
C
C         THE VALUE OF IFLAG SHOULD NOT BE CHANGED BY FCN UNLESS
C         THE USER WANTS TO TERMINATE EXECUTION OF HYBRD1.
C         IN THIS CASE SET IFLAG TO A NEGATIVE INTEGER.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF FUNCTIONS AND VARIABLES.
C
C       X IS AN ARRAY OF LENGTH N. ON INPUT X MUST CONTAIN
C         AN INITIAL ESTIMATE OF THE SOLUTION VECTOR. ON OUTPUT X
C         CONTAINS THE FINAL ESTIMATE OF THE SOLUTION VECTOR.
C
C       FVEC IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS
C         THE FUNCTIONS EVALUATED AT THE OUTPUT X.
C
C       TOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
C         WHEN THE ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
C         BETWEEN X AND THE SOLUTION IS AT MOST TOL.
C
C       INFO IS AN INTEGER OUTPUT VARIABLE. IF THE USER HAS
C         TERMINATED EXECUTION, INFO IS SET TO THE (NEGATIVE)
C         VALUE OF IFLAG. SEE DESCRIPTION OF FCN. OTHERWISE,
C         INFO IS SET AS FOLLOWS.
C
C         INFO = 0   IMPROPER INPUT PARAMETERS.
C
C         INFO = 1   ALGORITHM ESTIMATES THAT THE RELATIVE ERROR
C                    BETWEEN X AND THE SOLUTION IS AT MOST TOL.
C
C         INFO = 2   NUMBER OF CALLS TO FCN HAS REACHED OR EXCEEDED
C                    200*(N+1).
C
C         INFO = 3   TOL IS TOO SMALL. NO FURTHER IMPROVEMENT IN
C                    THE APPROXIMATE SOLUTION X IS POSSIBLE.
C
C         INFO = 4   ITERATION IS NOT MAKING GOOD PROGRESS.
C
C       WA IS A WORK ARRAY OF LENGTH LWA.
C
C       LWA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
C         (N*(3*N+13))/2.
C
C     SUBPROGRAMS CALLED
C
C       USER-SUPPLIED ...... FCN
C
C       MINPACK-SUPPLIED ... HYBRD
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER INDEX,J,LR,MAXFEV,ML,MODE,MU,NFEV,NPRINT
      DOUBLE PRECISION EPSFCN,FACTOR,ONE,XTOL,ZERO
      DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/
      INFO = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. TOL .LT. ZERO .OR. LWA .LT. (N*(3*N + 13))/2)
     *   GO TO 20
C
C     CALL HYBRD.
C
      MAXFEV = 200*(N + 1)
      XTOL = TOL
      ML = N - 1
      MU = N - 1
      EPSFCN = ZERO
      MODE = 2
c     DO 10 J = 1, N
      DO 10 J = 1, LWA
         WA(J) = ONE
   10    CONTINUE
      NPRINT = 0
      LR = (N*(N + 1))/2
      INDEX = 6*N + LR
      CALL HYBRD(FCN,N,X,FVEC,XTOL,MAXFEV,ML,MU,EPSFCN,WA(1),MODE,
     *           FACTOR,NPRINT,INFO,NFEV,WA(INDEX+1),N,WA(6*N+1),LR,
     *           WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1))
      IF (INFO .EQ. 5) INFO = 4
   20 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE HYBRD1.
C
      END
      SUBROUTINE QFORM(M,N,Q,LDQ,WA)
      INCLUDE 'EMMIX.max'
      INTEGER M,N,LDQ
      DOUBLE PRECISION Q(LDQ,M),WA(M)
C     **********
C
C     SUBROUTINE QFORM
C
C     THIS SUBROUTINE PROCEEDS FROM THE COMPUTED QR FACTORIZATION OF
C     AN M BY N MATRIX A TO ACCUMULATE THE M BY M ORTHOGONAL MATRIX
C     Q FROM ITS FACTORED FORM.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE QFORM(M,N,Q,LDQ,WA)
C
C     WHERE
C
C       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF ROWS OF A AND THE ORDER OF Q.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF COLUMNS OF A.
C
C       Q IS AN M BY M ARRAY. ON INPUT THE FULL LOWER TRAPEZOID IN
C         THE FIRST MIN(M,N) COLUMNS OF Q CONTAINS THE FACTORED FORM.
C         ON OUTPUT Q HAS BEEN ACCUMULATED INTO A SQUARE MATRIX.
C
C       LDQ IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY Q.
C
C       WA IS A WORK ARRAY OF LENGTH M.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... MIN0
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,J,JM1,K,L,MINMN,NP1
      DOUBLE PRECISION ONE,SUM,TEMP,ZERO
      DATA ONE,ZERO /1.0D0,0.0D0/
C
C     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
C
      MINMN = MIN0(M,N)
      IF (MINMN .LT. 2) GO TO 30
      DO 20 J = 2, MINMN
         JM1 = J - 1
         DO 10 I = 1, JM1
            Q(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
C
      NP1 = N + 1
      IF (M .LT. NP1) GO TO 60
      DO 50 J = NP1, M
         DO 40 I = 1, M
            Q(I,J) = ZERO
   40       CONTINUE
         Q(J,J) = ONE
   50    CONTINUE
   60 CONTINUE
C
C     ACCUMULATE Q FROM ITS FACTORED FORM.
C
      DO 120 L = 1, MINMN
         K = MINMN - L + 1
         DO 70 I = K, M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   70       CONTINUE
         Q(K,K) = ONE
         IF (WA(K) .EQ. ZERO) GO TO 110
         DO 100 J = K, M
            SUM = ZERO
            DO 80 I = K, M
               SUM = SUM + Q(I,J)*WA(I)
   80          CONTINUE
            TEMP = SUM/WA(K)
            DO 90 I = K, M
               Q(I,J) = Q(I,J) - TEMP*WA(I)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE QFORM.
C
      END
      SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
      INCLUDE 'EMMIX.max'
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(LIPVT)
      LOGICAL PIVOT
      DOUBLE PRECISION A(LDA,MAXNG),RDIAG(MAXNG),ACNORM(MAXNG),WA(MAXNG)
C     **********
C
C     SUBROUTINE QRFAC
C
C     THIS SUBROUTINE USES HOUSEHOLDER TRANSFORMATIONS WITH COLUMN
C     PIVOTING (OPTIONAL) TO COMPUTE A QR FACTORIZATION OF THE
C     M BY N MATRIX A. THAT IS, QRFAC DETERMINES AN ORTHOGONAL
C     MATRIX Q, A PERMUTATION MATRIX P, AND AN UPPER TRAPEZOIDAL
C     MATRIX R WITH DIAGONAL ELEMENTS OF NONINCREASING MAGNITUDE,
C     SUCH THAT A*P = Q*R. THE HOUSEHOLDER TRANSFORMATION FOR
C     COLUMN K, K = 1,2,...,MIN(M,N), IS OF THE FORM
C
C                           T
C           I - (1/U(K))*U*U
C
C     WHERE U HAS ZEROS IN THE FIRST K-1 POSITIONS. THE FORM OF
C     THIS TRANSFORMATION AND THE METHOD OF PIVOTING FIRST
C     APPEARED IN THE CORRESPONDING LINPACK SUBROUTINE.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,RDIAG,ACNORM,WA)
C
C     WHERE
C
C       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF ROWS OF A.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF COLUMNS OF A.
C
C       A IS AN M BY N ARRAY. ON INPUT A CONTAINS THE MATRIX FOR
C         WHICH THE QR FACTORIZATION IS TO BE COMPUTED. ON OUTPUT
C         THE STRICT UPPER TRAPEZOIDAL PART OF A CONTAINS THE STRICT
C         UPPER TRAPEZOIDAL PART OF R, AND THE LOWER TRAPEZOIDAL
C         PART OF A CONTAINS A FACTORED FORM OF Q (THE NON-TRIVIAL
C         ELEMENTS OF THE U VECTORS DESCRIBED ABOVE).
C
C       LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.
C
C       PIVOT IS A LOGICAL INPUT VARIABLE. IF PIVOT IS SET TRUE,
C         THEN COLUMN PIVOTING IS ENFORCED. IF PIVOT IS SET FALSE,
C         THEN NO COLUMN PIVOTING IS DONE.
C
C       IPVT IS AN INTEGER OUTPUT ARRAY OF LENGTH LIPVT. IPVT
C         DEFINES THE PERMUTATION MATRIX P SUCH THAT A*P = Q*R.
C         COLUMN J OF P IS COLUMN IPVT(J) OF THE IDENTITY MATRIX.
C         IF PIVOT IS FALSE, IPVT IS NOT REFERENCED.
C
C       LIPVT IS A POSITIVE INTEGER INPUT VARIABLE. IF PIVOT IS FALSE,
C         THEN LIPVT MAY BE AS SMALL AS 1. IF PIVOT IS TRUE, THEN
C         LIPVT MUST BE AT LEAST N.
C
C       RDIAG IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         DIAGONAL ELEMENTS OF R.
C
C       ACNORM IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         NORMS OF THE CORRESPONDING COLUMNS OF THE INPUT MATRIX A.
C         IF THIS INFORMATION IS NOT NEEDED, THEN ACNORM CAN COINCIDE
C         WITH RDIAG.
C
C       WA IS A WORK ARRAY OF LENGTH N. IF PIVOT IS FALSE, THEN WA
C         CAN COINCIDE WITH RDIAG.
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... DPMPAR,ENORM
C
C       FORTRAN-SUPPLIED ... DMAX1,DSQRT,MIN0
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,J,JP1,K,KMAX,MINMN
      DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      DOUBLE PRECISION DPMPAR,ENORM
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
      EPSMCH = DPMPAR(1)
C
C     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
C
      DO 10 J = 1, N
         ACNORM(J) = ENORM(M,A(1,J))
         RDIAG(J) = ACNORM(J)
         WA(J) = RDIAG(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
C
C     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
C
      MINMN = MIN0(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
C
C        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
C
         KMAX = J
         DO 20 K = J, N
            IF (RDIAG(K) .GT. RDIAG(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .EQ. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         RDIAG(KMAX) = RDIAG(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
C
C        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
C        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
C
         AJNORM = ENORM(M-J+1,A(J,J))
         IF (AJNORM .EQ. ZERO) GO TO 100
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
C
C        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
C        AND UPDATE THE NORMS.
C
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. RDIAG(K) .EQ. ZERO) GO TO 80
            TEMP = A(J,K)/RDIAG(K)
            RDIAG(K) = RDIAG(K)*DSQRT(DMAX1(ZERO,ONE-TEMP**2))
            IF (P05*(RDIAG(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
            RDIAG(K) = ENORM(M-J,A(JP1,K))
            WA(K) = RDIAG(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         RDIAG(J) = -AJNORM
  110    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE QRFAC.
C
      END
      SUBROUTINE R1MPYQ(M,N,A,LDA,V,W)
      INCLUDE 'EMMIX.max'
      INTEGER M,N,LDA
      DOUBLE PRECISION A(LDA,MAXNG),V(MAXNG),W(MAXNG)
C     **********
C
C     SUBROUTINE R1MPYQ
C
C     GIVEN AN M BY N MATRIX A, THIS SUBROUTINE COMPUTES A*Q WHERE
C     Q IS THE PRODUCT OF 2*(N - 1) TRANSFORMATIONS
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     AND GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE WHICH
C     ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES, RESPECTIVELY.
C     Q ITSELF IS NOT GIVEN, RATHER THE INFORMATION TO RECOVER THE
C     GV, GW ROTATIONS IS SUPPLIED.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE R1MPYQ(M,N,A,LDA,V,W)
C
C     WHERE
C
C       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF ROWS OF A.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF COLUMNS OF A.
C
C       A IS AN M BY N ARRAY. ON INPUT A MUST CONTAIN THE MATRIX
C         TO BE POSTMULTIPLIED BY THE ORTHOGONAL MATRIX Q
C         DESCRIBED ABOVE. ON OUTPUT A*Q HAS REPLACED A.
C
C       LDA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN M
C         WHICH SPECIFIES THE LEADING DIMENSION OF THE ARRAY A.
C
C       V IS AN INPUT ARRAY OF LENGTH N. V(I) MUST CONTAIN THE
C         INFORMATION NECESSARY TO RECOVER THE GIVENS ROTATION GV(I)
C         DESCRIBED ABOVE.
C
C       W IS AN INPUT ARRAY OF LENGTH N. W(I) MUST CONTAIN THE
C         INFORMATION NECESSARY TO RECOVER THE GIVENS ROTATION GW(I)
C         DESCRIBED ABOVE.
C
C     SUBROUTINES CALLED
C
C       FORTRAN-SUPPLIED ... DABS,DSQRT
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,J,NMJ,NM1
      DOUBLE PRECISION COS,ONE,SIN,TEMP
      DATA ONE /1.0D0/
C
C     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 50
      DO 20 NMJ = 1, NM1
         J = N - NMJ
         IF (DABS(V(J)) .GT. ONE) COS = ONE/V(J)
         IF (DABS(V(J)) .GT. ONE) SIN = DSQRT(ONE-COS**2)
         IF (DABS(V(J)) .LE. ONE) SIN = V(J)
         IF (DABS(V(J)) .LE. ONE) COS = DSQRT(ONE-SIN**2)
         DO 10 I = 1, M
            TEMP = COS*A(I,J) - SIN*A(I,N)
            A(I,N) = SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
C
C     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
C
      DO 40 J = 1, NM1
         IF (DABS(W(J)) .GT. ONE) COS = ONE/W(J)
         IF (DABS(W(J)) .GT. ONE) SIN = DSQRT(ONE-COS**2)
         IF (DABS(W(J)) .LE. ONE) SIN = W(J)
         IF (DABS(W(J)) .LE. ONE) COS = DSQRT(ONE-SIN**2)
         DO 30 I = 1, M
            TEMP = COS*A(I,J) + SIN*A(I,N)
            A(I,N) = (-SIN)*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE R1MPYQ.
C
      END
      SUBROUTINE R1UPDT(M,N,S,LS,U,V,W,SING)
      INCLUDE 'EMMIX.max'
      INTEGER M,N,LS
      LOGICAL SING
      DOUBLE PRECISION S(LS),U(MAXNG),V(MAXNG),W(MAXNG)
C     **********
C
C     SUBROUTINE R1UPDT
C
C     GIVEN AN M BY N LOWER TRAPEZOIDAL MATRIX S, AN M-VECTOR U,
C     AND AN N-VECTOR V, THE PROBLEM IS TO DETERMINE AN
C     ORTHOGONAL MATRIX Q SUCH THAT
C
C                   T
C           (S + U*V )*Q
C
C     IS AGAIN LOWER TRAPEZOIDAL.
C
C     THIS SUBROUTINE DETERMINES Q AS THE PRODUCT OF 2*(N - 1)
C     TRANSFORMATIONS
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     WHERE GV(I), GW(I) ARE GIVENS ROTATIONS IN THE (I,N) PLANE
C     WHICH ELIMINATE ELEMENTS IN THE I-TH AND N-TH PLANES,
C     RESPECTIVELY. Q ITSELF IS NOT ACCUMULATED, RATHER THE
C     INFORMATION TO RECOVER THE GV, GW ROTATIONS IS RETURNED.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE R1UPDT(M,N,S,LS,U,V,W,SING)
C
C     WHERE
C
C       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF ROWS OF S.
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
C         OF COLUMNS OF S. N MUST NOT EXCEED M.
C
C       S IS AN ARRAY OF LENGTH LS. ON INPUT S MUST CONTAIN THE LOWER
C         TRAPEZOIDAL MATRIX S STORED BY COLUMNS. ON OUTPUT S CONTAINS
C         THE LOWER TRAPEZOIDAL MATRIX PRODUCED AS DESCRIBED ABOVE.
C
C       LS IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
C         (N*(2*M-N+1))/2.
C
C       U IS AN INPUT ARRAY OF LENGTH M WHICH MUST CONTAIN THE
C         VECTOR U.
C
C       V IS AN ARRAY OF LENGTH N. ON INPUT V MUST CONTAIN THE VECTOR
C         V. ON OUTPUT V(I) CONTAINS THE INFORMATION NECESSARY TO
C         RECOVER THE GIVENS ROTATION GV(I) DESCRIBED ABOVE.
C
C       W IS AN OUTPUT ARRAY OF LENGTH M. W(I) CONTAINS INFORMATION
C         NECESSARY TO RECOVER THE GIVENS ROTATION GW(I) DESCRIBED
C         ABOVE.
C
C       SING IS A LOGICAL OUTPUT VARIABLE. SING IS SET TRUE IF ANY
C         OF THE DIAGONAL ELEMENTS OF THE OUTPUT S ARE ZERO. OTHERWISE
C         SING IS SET FALSE.
C
C     SUBPROGRAMS CALLED
C
C       MINPACK-SUPPLIED ... DPMPAR
C
C       FORTRAN-SUPPLIED ... DABS,DSQRT
C
C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE,
C     JOHN L. NAZARETH
C
C     **********
      INTEGER I,J,JJ,L,NMJ,NM1
      DOUBLE PRECISION COS,COTAN,GIANT,ONE,P5,P25,SIN,TAN,TAU,TEMP,
     *                 ZERO
      DOUBLE PRECISION DPMPAR
      DATA ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/
C
C     GIANT IS THE LARGEST MAGNITUDE.
C
      GIANT = DPMPAR(3)
C
C     INITIALIZE THE DIAGONAL ELEMENT POINTER.
C
      JJ = (N*(2*M - N + 1))/2 - (M - N)
C
C     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.
C
      L = JJ
      DO 10 I = N, M
         W(I) = S(L)
         L = L + 1
   10    CONTINUE
C
C     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
C     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 NMJ = 1, NM1
         J = N - NMJ
         JJ = JJ - (M - J + 1)
         W(J) = ZERO
         IF (V(J) .EQ. ZERO) GO TO 50
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF V.
C
         IF (DABS(V(N)) .GE. DABS(V(J))) GO TO 20
            COTAN = V(N)/V(J)
            SIN = P5/DSQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (DABS(COS)*GIANT .GT. ONE) TAU = ONE/COS
            GO TO 30
   20    CONTINUE
            TAN = V(J)/V(N)
            COS = P5/DSQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
   30    CONTINUE
C
C        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
C        NECESSARY TO RECOVER THE GIVENS ROTATION.
C
         V(N) = SIN*V(J) + COS*V(N)
         V(J) = TAU
C
C        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
C
         L = JJ
         DO 40 I = J, M
            TEMP = COS*S(L) - SIN*W(I)
            W(I) = SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.
C
      DO 80 I = 1, M
         W(I) = W(I) + V(N)*U(I)
   80    CONTINUE
C
C     ELIMINATE THE SPIKE.
C
      SING = .FALSE.
      IF (NM1 .LT. 1) GO TO 140
      DO 130 J = 1, NM1
         IF (W(J) .EQ. ZERO) GO TO 120
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF THE SPIKE.
C
         IF (DABS(S(JJ)) .GE. DABS(W(J))) GO TO 90
            COTAN = S(JJ)/W(J)
            SIN = P5/DSQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (DABS(COS)*GIANT .GT. ONE) TAU = ONE/COS
            GO TO 100
   90    CONTINUE
            TAN = W(J)/S(JJ)
            COS = P5/DSQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
  100    CONTINUE
C
C        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
C
         L = JJ
         DO 110 I = J, M
            TEMP = COS*S(L) + SIN*W(I)
            W(I) = (-SIN)*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
  110       CONTINUE
C
C        STORE THE INFORMATION NECESSARY TO RECOVER THE
C        GIVENS ROTATION.
C
         W(J) = TAU
  120    CONTINUE
C
C        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
C
         IF (S(JJ) .EQ. ZERO) SING = .TRUE.
         JJ = JJ + (M - J + 1)
  130    CONTINUE
  140 CONTINUE
C
C     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.
C
      L = JJ
      DO 150 I = N, M
         S(L) = W(I)
         L = L + 1
  150    CONTINUE
      IF (S(JJ) .EQ. ZERO) SING = .TRUE.
      RETURN
C
C     LAST CARD OF SUBROUTINE R1UPDT.
C
      END

C
C
C  This group of subroutines implements the K-means Clustering algorithm.
C  Implemented by David Peel May 1994                      

      SUBROUTINE KMEANS(NIND,NATT,NG,X,IDT,EPSILON,TT,IER)
C      Main subroutine
      
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INTEGER T,TT
       EXTERNAL RANDNUM
       INTEGER FLAGS(40),FYLENO
       COMMON /STORE2/ FLAGS,FYLENO
       DOUBLE PRECISION RANDNUM
       INCLUDE 'EMMIX.max'
       DIMENSION X(MNIND,MNATT),XK(MAXNG,MNATT),
     &           XKOLD(MAXNG,MNATT),IDT(MNIND),
     &           XSTAN(MNIND,MNATT)
      IER=0
      CALL KSTAND(NIND,NATT,X,XSTAN) 
      CALL KSEED(NIND,NATT,NG,XSTAN,XK,IER)
      DO 30 T=1,TT
      DO 430 KK=1,NG
        DO 420 LL=1,NATT 
            XKOLD(KK,LL)=XK(KK,LL)
420   	  CONTINUE
430     CONTINUE
	DO 20 KK=1,NIND
          CALL WINNER(NATT,NG,XSTAN,KK,XK,IDT,IER)
20      CONTINUE
          CALL UPDATE(NIND,NATT,NG,XSTAN,XK,IDT,IER)
          ET=RULE(NG,NATT,XKOLD,XK)
          IF (ET.LE.EPSILON) GO TO 99
30    CONTINUE
      
      WRITE (FYLENO,*) 'REACHED MAXIMUM NUMBER OF ',TT,' ITERATIONS'
      IER=-41
99    RETURN
      END

      SUBROUTINE KSTAND(NIND,NATT,X,XNEW) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      DIMENSION X(MNIND,MNATT),XNEW(MNIND,MNATT),
     &          XVAR(MNATT),XMU(MNATT)
      DO 200 J=1,NATT
       XMU(J)=0
       DO 200 I=1,NIND
        XMU(J)=XMU(J)+X(I,J)/NIND 
200   CONTINUE       
      DO 210 J=1,NATT
       XVAR(J)=0
       DO 210 I=1,NIND
        XVAR(J)=XVAR(J)+(X(I,J)-XMU(J))*
     &            (X(I,J)-XMU(J))/(NIND-1)
210   CONTINUE
      DO 220 J=1,NATT
       DO 220 I=1,NIND
         XNEW(I,J)=(X(I,J)-XMU(J))/XVAR(J)
220   CONTINUE
      RETURN
      END

      SUBROUTINE KSEED(NIND,NATT,NG,XSTAN,XK,IER)               
c     This Subroutine chooses the initial K seeds (Means of clusters)
c     for the algorithm. At present they are chosen from data set at 
c     random.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER CHOICE
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM
       INCLUDE 'EMMIX.max'
      DIMENSION XSTAN(MNIND,MNATT),XK(MAXNG,MNATT)
      DO 210 I=1,NG
        R=RANDNUM()
        R=R*NIND
c       Convert CHOICE to integer
	CHOICE=INT(R)+1
	DO 200 J=1,NATT
          XK(I,J)=XSTAN(CHOICE,J) 
200	CONTINUE
210   CONTINUE
      RETURN
      END

      SUBROUTINE WINNER(NATT,NG,XSTAN,KK,XK,IDT,IER)
c     This subroutine determines the allocation of the KKth point 
c     ie which mean is closest to the given data point (Euclidean).

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      DIMENSION XSTAN(MNIND,MNATT),XK(MAXNG,MNATT),
     &          IDT(MNIND)
      DO 310 I=1,NG
        DIST=0
  	DO 300 J=1,NATT
          DIST=DIST+(XSTAN(KK,J)-XK(I,J))**2
300 	CONTINUE
        IF (I.EQ.1) DISTB=DIST 
        IF (DIST.LE.DISTB) THEN
	  IDT(KK)=I
          DISTB=DIST
        ENDIF
310   CONTINUE
      RETURN
      END 
          
      SUBROUTINE UPDATE(NIND,NATT,NG,XSTAN,XK,IDT,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      DIMENSION XK(MAXNG,MNATT),IDT(MNIND),
     &          XSTAN(MNIND,MNATT),N(MAXNG)
        DO 410 II=1,NG
          N(II)=0
          DO 410 LL=1,NATT
            XK(II,LL)=0 
410     CONTINUE
        DO 450 I=1,NIND
          II=IDT(I)
          N(II)=N(II)+1
c         Update rules
          DO 440 LL=1,NATT
            XK(II,LL)=XK(II,LL)+XSTAN(I,LL)
440       CONTINUE
450     CONTINUE
        DO 499 II=1,NG
          DO 499 LL=1,NATT
           IF (N(II).NE.0) THEN
            XK(II,LL)=XK(II,LL)/N(II)       
           ELSE
             RETURN
           ENDIF
499     CONTINUE
      RETURN
      END


      FUNCTION RULE(NG,NATT,XKOLD,XK)
c     This function returns the value used to determine if the algorithm
c     has converged it is a measure of the change in the nodes from iteration 
c     to iteration. 
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      INTEGER R
      DIMENSION XK(MAXNG,MNATT),XKOLD(MAXNG,MNATT)
      RULE=0.0
      DO 510 KK=1,NATT
        DO 500 R=1,NG
           RULE=RULE+abs(XK(R,KK)-XKOLD(R,KK))
500     CONTINUE
510   CONTINUE
      RETURN
      END

C
C
C
Copyright (C) 1985 Numerical Recipes Software -- GAMMLN
      FUNCTION GAMMLN(XX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION COF(6),STP,HALF,ONE,FPF,X,TMP,SER,XX
c      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
c    *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA COF,STP/76.18009173,-86.50532033,24.01409822,
     *    -1.231739516,.120858003E-2,-.536382E-5,2.50662827465/
      DATA HALF,ONE,FPF/0.5,1.0,5.5/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END

      FUNCTION DIGAMA(X, IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     ALGORITHM AS 103  APPL. STATIST. (1976) VOL.25, NO.3
C
C     Calculates DIGAMMA(X) = D( LOG( GAMMA(X))) / DX
C
      DOUBLE PRECISION  ZERO, HALF, ONE
C
C     Set constants, SN = Nth Stirling coefficient, D1 = DIGAMMA(1.0)
C
      DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/
      DATA S, C, S3, S4, S5, D1 /1.E-05, 8.5, 8.333333333E-02,
     *    8.3333333333E-03, 3.96825 3968E-03, -0.57721 56649/
C
C     Check argument is positive
C
      DIGAMA = ZERO
      Y = X
      IER = 1
      IF (Y .LE. ZERO) RETURN
      IER = 0
C
C     Use approximation if argument <= S
C
      IF (Y .LE. S) THEN
        DIGAMA = D1 - ONE / Y
        RETURN
      END IF
C
C     Reduce to DIGAMA(X + N) where (X + N) >= C
C
    1 IF (Y .GE. C) GO TO 2
      DIGAMA = DIGAMA - ONE/Y
      Y = Y + ONE
      GO TO 1
C
C     Use Stirling's (actually de Moivre's) expansion if argument > C
C
    2 R = ONE / Y
      DIGAMA = DIGAMA + LOG(Y) - HALF*R
      R = R * R
      DIGAMA = DIGAMA - R*(S3 - R*(S4 - R*S5))
      RETURN
      END

      SUBROUTINE TMOM(NIND,NATT,NG,X,XMU,XVAR,W,XUU)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      DIMENSION X(MNIND,MNATT),XMU(MAXNG,MNATT),XUU(MAXNG)
      DIMENSION XVAR(MAXNG,MNATT,MNATT),W(MNIND,MAXNG)
      DIMENSION C(2,MAXNG)
      DO 150 K=1,NG
      DO 100 I=1,NATT
        C(1,K)=C(1,K)+XVAR(K,I,I)*XVAR(K,I,I)
100   CONTINUE 
      CBOT=0
      C(2,K)=0
      DO 110 I=1,NIND
        DO 110 J=1,NATT
          C(2,K)=C(2,K)+W(I,K)*(X(I,J)-XMU(K,J))**4.0
          CBOT=CBOT+W(I,K)
110   CONTINUE
      C(2,K)=C(2,K)/CBOT
      XUU(K)=(6.0*C(1,K)-4.0*C(2,K))/(3*C(1,K)-C(2,K))
150   CONTINUE
      RETURN
      END

      SUBROUTINE TEQ0(NDUMMY,XUU1,FVEC1,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      COMMON /HYBD/ UTEMP(MNIND,MAXNG),TNIND,TNATT,XUUOLD(MAXNG),
     &   WTEMP(MNIND,MAXNG),NGTEMP
        FVEC1=0
      DO 110 K=1,NGTEMP
        DO 100 I=1,TNIND
          IF (XUU1.LE.0.00001)THEN
            XUU1=0.00001
            IER=1
            RETURN
           ENDIF
          IER=0
          DIG1=DIGAMA(0.5*XUU1,IER)
          IF (IER.GT.0) RETURN
          DIG2=DIGAMA(0.5*(XUUOLD(1)+TNATT),IER)
          IF (IER.GT.0) RETURN
          XLOG1=LOG(0.5*XUU1)
          XLOG2=LOG(0.5*(XUUOLD(1)+TNATT))
          XSUM=(LOG(UTEMP(I,K))-UTEMP(I,K))
          FVEC1=FVEC1+WTEMP(I,K)*((-1)*DIG1+XLOG1+1+XSUM+DIG2-XLOG2)
100     CONTINUE
110   CONTINUE
      RETURN
      END


      SUBROUTINE TEQ(NG,XUU,FVEC,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'EMMIX.max'
      COMMON /HYBD/ UTEMP(MNIND,MAXNG),TNIND,TNATT,XUUOLD(MAXNG),
     &   WTEMP(MNIND,MAXNG),NGTEMP
      DIMENSION XUU(MAXNG),FVEC(MAXNG)
      DO 110 K=1,NG 
        FVEC(K)=0
        DO 100 I=1,TNIND
          IF (XUU(K).LE.0.00001) THEN
          XUU(K)=0.00001
          IER=1
          RETURN
          ENDIF
          IER=0
          DIG1=DIGAMA(0.5*XUU(K),IER)
          IF (IER.GT.0) RETURN
          DIG2=DIGAMA(0.5*(XUUOLD(K)+TNATT),IER)
          IF (IER.GT.0) RETURN
          XLOG1=LOG(0.5*XUU(K))
          XLOG2=LOG(0.5*(XUUOLD(K)+TNATT))
          XSUM=(LOG(UTEMP(I,K))-UTEMP(I,K))
          FVEC(K)=FVEC(K)+WTEMP(I,K)*((-1)*DIG1+XLOG1+1+XSUM+DIG2-XLOG2)
100     CONTINUE
110   CONTINUE
      RETURN
      END

        SUBROUTINE TFREE(NIND,NATT,NG,XUU,U,W,ITER,IER)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INCLUDE 'EMMIX.max'
      COMMON /HYBD/ UTEMP(MNIND,MAXNG),TNIND,TNATT,XUUOLD(MAXNG),
     &    WTEMP(MNIND,MAXNG),NGTEMP
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO 
        EXTERNAL TEQ,TEQ0
c	common/comm6/fail,s0jeq1,nboot,group,itct
	DIMENSION fvec(MAXNG),W(MNIND,MAXNG)
c	integer fail,s0jeq1,nboot,group,itct
c	common/comm7/est,newest
c	double precision newest(8)
c	double precision est(8)

c	double precision beta,gamma
c	double precision newbetai,newgammai
	double precision xtol
	double precision work((MAXNG*(3*MAXNG+13))/2)
	integer sizehy,ifail

        DIMENSION XUU(MAXNG),U(MNIND,MAXNG)
        NGTEMP=NG
        DO 100 K=1,NG
         XUUOLD(K)=XUU(K)
         DO 100 I=1,NIND
          UTEMP(I,K)=U(I,K)
          WTEMP(I,K)=W(I,K)
100     CONTINUE           
        TNIND=FLOAT(NIND)
        TNATT=FLOAT(NATT)

C       	if(s0jeq1 .eq. 2) then
C	  indlo=2
C	else
C	  indlo=1
C	endif
C	do 20 group=indlo,2
C	  tmp=(group-1)*3+1
	  ifail=0
	  xtol=0.000000001
	  xtol=0.000000000001
	  sizehy=(NG*(3*NG+13))/2
          IF (FLAGS(7).EQ.3) THEN
          XUU1=XUU(1)
	  call HYBRD1(TEQ0,1,XUU1,FVEC1,xtol,ifail,work,sizehy)
      DO 111 K=1,NG 
       XUU(K)=XUU1
111   CONTINUE
          ELSE
	  call HYBRD1(TEQ,NG,XUU,FVEC,xtol,ifail,work,sizehy)
          ENDIF
C20     continue
      RETURN
      END
C
C
C   This group of subroutines deal with generating random numbers
C   
C   The subroutines in this file  basically act as a interface between
C   how this program calls the random number generator ( R=RANDNUM(SEED)
C  

C   D.Peel Nov 1995

       SUBROUTINE DETERRANDOM(IER)
C      This subroutine is called at the very beginning of the program
C      to determine if the random number generator is working 
C
C      The result of the test is stored in the common random 
C      variable RANDTYPE
C      This is now defunct as only one random number generator is tried
C      but is left in the code in case a number of generators are available
C      in future versions

	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	   COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
	   DIMENSION XISEEDS(1000)
           IER=0
	   DO 100 K=1,1000
             XISEEDS(K)=0
100    CONTINUE
           IX=10023
           IY=324
           IZ=54367
           XISEEDS(1)=.435543
           DO 106 I=2,1000
             XISEEDS(I)=RANDOM(IX,IY,IZ)
             IF (XISEEDS(I).EQ.0) GO TO 107
             DO 105 J=1,I-1
               IF (XISEEDS(I).EQ.XISEEDS(J)) GO TO 107
105      CONTINUE
106    CONTINUE
       RANDTYPE=1
           RETURN

107    WRITE(*,*)'PROBLEM: The random number generator'
       WRITE(*,*)'  does not seem to produce random numbers.'
       WRITE(*,*)'  Applied statistics random number generator failed'
       WRITE(*,*)'  to produce suitable random numbers as well. Modify'
       WRITE(*,*)'  the file nmmrand.f to make this program'
       WRITE(*,*)'  compatible with your inbuilt random number'
       WRITE(*,*)'  generator or alternatively use another random'
       WRITE(*,*)'  number generator.'
       WRITE(*,*)      
       WRITE(*,*)' At present MIXCLUS will still function but you will'
       WRITE(*,*)' be unable to use features that incorporate random'
       WRITE(*,*)' numbers'
       RANDTYPE=0
       IER=40
	   RETURN
	   END
	   
	   FUNCTION RANDNUM()
C          This is the function called by the program NMM. If you
C          wish to use your own portable random number generator 
C          then it should be used in place of this function.

	   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	   COMMON /STORE1/ SEED,RANDTYPE,IX,IY,IZ
	   IF (RANDTYPE.EQ.1) THEN
             RANDNUM=RANDOM(IX,IY,IZ)
	   ELSE
      WRITE(*,*)'ERROR: As previously described due to random number'
      WRITE(*,*)'       generator problems features utilising'
      WRITE(*,*)'       random numbers are unavailable'
	   STOP
	   ENDIF
           RETURN
	   END

      DOUBLE PRECISION FUNCTION RANDOM(IX,IY,IZ)
C
C     Algorithm AS 183 Appl. Statist. (1982) vol.31, no.2
C
C     Returns a pseudo-random number rectangularly distributed
C     between 0 and 1.   The cycle length is 6.95E+12 (See page 123
C     of Applied Statistics (1984) vol.33), not as claimed in the
C     original article.
C
C     IX, IY and IZ should be set to integer values between 1 and
C     30000 before the first entry.
C
C     Integer arithmetic up to 30323 is required.
C
      INTEGER IX, IY, IZ
c      COMMON /RANDC/ IX, IY, IZ

C
      IX = 171 * MOD(IX, 177) - 2 * (IX / 177)
      IY = 172 * MOD(IY, 176) - 35 * (IY / 176)
      IZ = 170 * MOD(IZ, 178) - 63 * (IZ / 178)
C
      IF (IX .LT. 0) IX = IX + 30269
      IF (IY .LT. 0) IY = IY + 30307
      IF (IZ .LT. 0) IZ = IZ + 30323
C
C     If integer arithmetic up to 5212632 is available, the preceding
C     6 statements may be replaced by:
C
C     IX = MOD(171 * IX, 30269)
C     IY = MOD(172 * IY, 30307)
C     IZ = MOD(170 * IZ, 30323)
C
      RANDOM = MOD(FLOAT(IX) / 30269. + FLOAT(IY) / 30307. +
     +                        FLOAT(IZ) / 30323., 1.0)
      RETURN
      END
C
C
C
C     The file contain subroutines contained in this file generate
C     bootstrap samples via three methods:
C      1/ Sampling from the original data set with replacement which
C         corresponds to multinomial weights on the data.
C      2/ A weighted likelihood approach which corresponds to
C         exponential weights on the data
C      3/ Parametric sample which is simply generating a sample from the
C         fitted estimates of the model on the original data set(normal
C         mixture model) 
C
      SUBROUTINE BSAMP(NIND,NATT,NG,XMU,XVAR,T,TX,X,WL,IDT,BT,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM
C
C     see head of driver file regarding this line
       INCLUDE 'EMMIX.max'
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      PARAMETER (MNIND=900)
C         maximum number of data points
c      PARAMETER (MNATT=20)
C         maximum dimensionality of data points
c      PARAMETER (MAXNG=20)
C         maximum number of groups
C          section.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION XVAR(MAXNG,MNATT,MNATT),
     &          T(MAXNG),TX(MNIND,MNATT),WL(MNIND),
     &          X(MNIND,MNATT),IDT(MNIND),XMU(MAXNG,MNATT)
      IER=0
      DO 10 I=1,NIND
        WL(I)=0.0
10    CONTINUE      
      IF ((FLAGS(18).EQ.1).AND.(BT.EQ.0)) THEN
       CALL SWR(NIND,WL)
      ELSEIF ((FLAGS(18).EQ.2).AND.(BT.EQ.0)) THEN 
       CALL WLIKE(NIND,WL)
      ELSE
       CALL MVNG(NIND,NATT,NG,XMU,XVAR,T,X,IDT,IER)
       DO 50 I=1,NIND
        WL(I)=1
50     CONTINUE
      ENDIF
      IF (FLAGS(11).EQ.2) THEN
        WRITE (FYLENO,*) '*************************************'
        WRITE (FYLENO,*) '* The Bootstrap Sample has NOT been *'
        WRITE (FYLENO,*) '* written since the space efficient *'
        WRITE (FYLENO,*) '* option has been chosen           *'
        WRITE (FYLENO,*) '*************************************'
      ENDIF
      RETURN
      END

      SUBROUTINE SWR(NIND,WL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RANDNUM
C
C     see head of driver file regarding this line
       INCLUDE 'EMMIX.max'
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      PARAMETER (MNIND=900)
C         maximum number of data points
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FLAGS(40),FYLENO,CHOICE
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION WL(MNIND)
      DOUBLE PRECISION RANDNUM
      DO 100 I=1,NIND
        R=RANDNUM()
        R=R*NIND
c       Convert CHOICE to integer
        CHOICE=INT(R)+1
          WL(CHOICE)=WL(CHOICE)+1
100   CONTINUE

      IF (FLAGS(11).NE.2) THEN
      WRITE (FYLENO,*) '   Bootstrap sample from last replicate' 
      WRITE (FYLENO,*) '     using sampling with replacement'
      WRITE (FYLENO,*) '    Index   Weight' 
      DO 110 I=1,NIND
        WRITE (FYLENO,*) I,WL(I)
110   CONTINUE
      ENDIF
      RETURN
      END
       

      SUBROUTINE WLIKE(NIND,WL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RANDNUM
C
C     see head of driver file regarding this line
       INCLUDE 'EMMIX.max'
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      PARAMETER (MNIND=900)
C         maximum number of data points
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION WL(MNIND)
      DOUBLE PRECISION RANDNUM
      WLSUM=0
      DO 100 I=1,NIND
        WEIGHT=RANDNUM()
c       Convert CHOICE to integer
        WEIGHT=log(WEIGHT) 
          WL(I)=WEIGHT
          WLSUM=WLSUM+WL(I)
100   CONTINUE
      DO 110 I=1,NIND
        WL(I)=WL(I)*NIND/WLSUM
110   CONTINUE
      IF (FLAGS(11).NE.2) THEN
      WRITE (FYLENO,*) '   Bootstrap sample from last replicate' 
      WRITE (FYLENO,*) '     using weighted likelihood approach' 
      WRITE (FYLENO,*) '    Index   Weight' 
      DO 120 I=1,NIND
        WRITE (FYLENO,*) I,WL(I)
120   CONTINUE
      ENDIF
      RETURN
      END



C     The group subroutines contained in this file generate random samples
C     from a multivariate normal mixture model
C     Written by K. Basford, modified by M. Latter, 1987. 
C     Further modified by Peter Adams Feb 1992 
C     for output to be written to pmvng.out1, .out2 and .out3 
C     Furthermore, nind, natt and ngp are NO LONGER written to pmvng.out2.
C     Modified by D.Peel SEP 1994 minor modifications and converted to a
C     subroutine to interact with MMresamp

      SUBROUTINE MVNG(NIND,NATT,NG,XMU,XVAR,T,X,IDT,IER)
C     This is the main subroutine to generate samples from a mixture of 
C     multivariate normals. The program uses the Naylor version to 
C     generate the multivariate normal data.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM
C
C     see head of driver file regarding this line
       INCLUDE 'EMMIX.max'
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      PARAMETER (MNIND=900)
C         maximum number of data points
c      PARAMETER (MNATT=20)
C         maximum dimensionality of data points
c      PARAMETER (MAXNG=20)
C         maximum number of groups
C          section.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION XVAR(MAXNG,MNATT,MNATT),XSAR(MNIND,MNATT),
     &          NT(MAXNG),T(MAXNG),
     &          X(MNIND,MNATT),IDT(MNIND),XMU(MAXNG,MNATT)
       IER=0
C      Determine the number in each group
        CALL RANNT(NIND,NG,T,NT,R,IER)
        IF (IER.NE.0) RETURN
C       Array NT contains the number of each group in the sample
        XKN=0     
        DO 30 K=1,NG
C         Generate all samples from each group at once
C         allocate to the y array in order
          KK=K
          CALL MVNORM(NT(K),NATT,KK,XMU,XVAR,XSAR,R)
C         Use group sample means and variances for initial allocation
          DO 20 I=1,NT(K)
          DO 20 J=1,NATT
            II=I+XKN
            X(II,J)=XSAR(I,J)
            IDT(II)=K
20        CONTINUE
          XKN=XKN+NT(K)
30      CONTINUE
32     FORMAT (2X,10F12.4)
33     FORMAT (10I4)
      IF (FLAGS(8).EQ.0) THEN
c       WRITE (FYLENO,*) '   Generated sample from MIXCLUS'
        DO 40 I=1,NIND
40        WRITE (FYLENO,32) (X(I,J),J=1,NATT)
        WRITE (FYLENO,33) (IDT(I),I=1,NIND)
      ELSEIF (FLAGS(11).NE.2) THEN
      WRITE (FYLENO,*) '   Bootstrap sample from last replicate' 
      WRITE (FYLENO,*) '     using parametric sampling' 
      WRITE (FYLENO,*) '       DATA-'

        DO 41 I=1,NIND
41        WRITE (FYLENO,32) (X(I,J),J=1,NATT)
       WRITE (FYLENO,*) '                  True Partition of Sample'
       WRITE (FYLENO,33) (IDT(I),I=1,NIND)
      ENDIF
99    CONTINUE

      RETURN 
      END


      SUBROUTINE RANNT(NIND,NG,T,NT,R,IER)
C     The samples have the number from each group
C     determined randomly in the same proportion  
C     as those assigned by the mixml program.
C     This subroutine written by D.Peel Dec 1994 differs from original
C     as no restriction on the number of groups. Which enables a simple
C     modification within harness to extend maximum Number of groups.

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     see head of driver file regarding this line
       INCLUDE 'EMMIX.max'
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      PARAMETER (MAXNG=20)
C         maximum number of groups
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER FLAGS(40),FYLENO
      COMMON /STORE2/ FLAGS,FYLENO
      DIMENSION T(MAXNG),NT(MAXNG),XL(MAXNG)
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM,R
      TEMP=0
      DO 110 K=1,NG
        NT(K)=0
        TEMP=TEMP+T(K)
        XL(K)=TEMP
110   CONTINUE
      DO 120 I=1,NIND
        R=RANDNUM()
        DO 115 K=1,NG
          IF (R.LE.XL(K)) THEN
            NT(K)=NT(K)+1
            GOTO 120
          ENDIF
115      CONTINUE
120    CONTINUE
       DO 130 I=1,NG
         IF (NT(I).LT.2) THEN
           IER=23
      WRITE(FYLENO,*) 'ERROR : Not enough points generated in group ',I
           RETURN
         ENDIF
130   CONTINUE
      RETURN
      END


      SUBROUTINE MVNORM(NS,NATT,KGP,XMU,XVAR,XSAR,R)
c      Multivariate normal sample generator

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     see head of driver file regarding this line
      INCLUDE 'EMMIX.max'
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      PARAMETER (MNIND=900)
C         maximum number of data points
c     PARAMETER (MNATT=20)
C         maximum dimensionality of data points
c     PARAMETER (MAXNG=20)
C         maximum number of groups
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DIMENSION X(MNATT),XVAR(MAXNG,MNATT,MNATT),
     &          C(MNATT,MNATT),XMU(MAXNG,MNATT),
     &          XSAR(MNIND,MNATT),Z(MNATT)

c      Initialize mean and variance arrays
       DO 11 I=1,NATT
         DO 11 J=1,NATT
          C(I,J)=0.0
11     CONTINUE
       DO 12 I=1,NATT
12      C(I,1)=XVAR(KGP,I,1)/SQRT(XVAR(KGP,1,1))
       DO 15 I=2,NATT
         IF(I.EQ.2) GO TO 177
         DO 13 J=2,I-1
           SUM=0.0
           DO 16 K=1,J-1
16           SUM=SUM+C(I,K)*C(J,K)    
           C(I,J)=(XVAR(KGP,I,J)-SUM)/C(J,J)
13       CONTINUE
177      CONTINUE
         SUM=0.0
         DO 14 K=1,I-1
14         SUM=SUM+C(I,K)*C(I,K)
         C(I,I)=SQRT(XVAR(KGP,I,I)-SUM)
15     CONTINUE
       DO 50 L=1,NS
         DO 17 I=1,NATT
          CALL NORM(Z(I),R)
17       CONTINUE
         DO 20 I=1,NATT
           SUM=0.0
           DO 18 J=1,I
18           SUM=SUM+C(I,J)*Z(J)
           X(I)=SUM+XMU(KGP,I)
20       CONTINUE
         DO 21 I=1,NATT
           XSAR(L,I)=X(I)
21       CONTINUE
50     CONTINUE
      RETURN 
      END
 



      SUBROUTINE NORM(Q,R)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL X,Y,G3X,V1,V2,SS,CL,Z
      EXTERNAL RANDNUM
      DOUBLE PRECISION RANDNUM,R
      R=RANDNUM()
      IF(R.LE..8638D0) THEN
10       R=RANDNUM()
         X=R
         R=RANDNUM()
         X=R+X
         R=RANDNUM()
         X=X+R
         Z=2.0D0*(X-1.50D0)
      ELSEIF(R.LE.0.9745D0) THEN
         R=RANDNUM()
         X=R
         R=RANDNUM()
         X=X+R
         Z=1.5D0*(X-1.D0)
      ELSEIF(R.LE.0.9973002039D0) THEN
20       R=RANDNUM()
         X=6.0D0*R-3.0D0
         R=X
         Y=0.358D0*RANDNUM()
         R=X
         IF(ABS(X).LT.1.0D0) THEN
          G3X=17.49731196D0*DEXP((-0.5D0)*X*X)-4.73570326D0*(3D0-X*X)
     &        -2.15787544D0*(1.5D0-ABS(X))
         ELSEIF(ABS(X).LT.1.5D0) THEN
      G3X=17.49731196D0*DEXP((-0.5D0)*X*X)-2.36785163D0*(3.0D0-ABS(X))
     &        **2.0D0-2.15787544D0*(1.5D0-ABS(X))
         ELSEIF(ABS(X).LT.3.0D0) THEN
      G3X=17.49731196D0*DEXP((-0.5D0)*X*X)-2.36785163D0*(3.0D0-ABS(X))
     &       **2.0D0
         ELSE
          G3X=0.0D0
         ENDIF
         IF(Y.LT.G3X)GO TO 30
         GO TO 20
30       CONTINUE
         Z=X
        ELSE
40      V1=RANDNUM()
        V2=RANDNUM()
        V1=2.0D0*V1-1.0D0
        V2=2.0D0*V2-1.0D0
        SS=V1*V1+V2*V2
        IF(SS.GT.1.0D0)GO TO 40
        CL=SQRT(9.0-2*LOG(SS)/SS)
        X=CL*V1
        Y=CL*V2
        IF(ABS(X).GT.3.0D0) THEN
          Z=X
          GOTO 45
        ELSEIF(ABS(Y).GT.3.0D0) THEN
          Z=Y
          Q=Z
          RETURN
        ENDIF
        GO TO 40
       ENDIF
45     Q=Z
      RETURN
      END
C
C
C
C Subroutine based on the program HACLUS by I.Delacy U.Q 1981                   
C Program to perform hierarchical agglomerative clustering.              
C implements the algorithm of Wishart 1969 Biometrics.                 
C clustering can be performed on:                                       
C     either variance standardized or unstandardised data,              
C Missing data is indicated by a value .LE.-90.0.                       
C therefore data must be scaled to .LE.ABS(89.9)                         
C                                                                       
C    The original program was written for general use as a clustering tool
C    with many extra options that are redundant when used in this context 
C    so they have been removed eg optional data input as a dissimilarity matrix.
C
C These subroutines will most probably be replaced in future versions.


      SUBROUTINE HIER(NIND,NATT,NG,X,IDT,ISU,IS,BETA,IFAULT)
C      The main controlling S/R for hierarchical clustering program          
C      Parameters:                                                           
C      IFIN : Y  finish; N  repeat cycle(reread data & cluster again)      
C      CUT   IGE  : G  Clusters genotypes; E  Clusters environments. 
C      ISU  : S  Standardize;  U  Don't standardize.                       
C      IS   : Clustering strategy (see S/R Disq).                          
C      BET  : BETA parameter for flexible sorting strategy.                
C      NIND   : NO. of Taxa. (Observations)                                    
C      NATT   : NO. of Attributes.                                           

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL CHECK
       INCLUDE 'EMMIX.max'
      DIMENSION Z(LIMZ),X(MNIND,MNATT),IDT(MNIND) 
      CHARACTER*5    H1(2,3),H2(2,3),H3(2,3),HH1(3),HH2(2)
      CHARACTER*5    HH4(5),HD(6),HS(7,6),FORM1(14),FORM3(11)
C      Set up character arrays storing strings naming and describing the 
C      various clustering methods used

      DATA (H1(1,J),J=1,3)/'GENOT','YPES ','     '/  
      DATA (H1(2,J),J=1,3)/'ENVIR','ONMEN','TS   '/ 
      DATA (H2(1,J),J=1,3)/'UNSTA','NDARD','IZED '/
      DATA (H2(2,J),J=1,3)/'STAND','ARDIZ','ED   '/
      DATA (H3(1,J),J=1,3)/'MEANS','     ','     '/
      DATA (H3(2,J),J=1,3)/'GXE E','FFECT','S    '/
      DATA HH1/'GROUP','S WIT','H    '/
      DATA HH2/'DATA ','ON   '/
      DATA HH4/'DISSI','MILAR','ITY M','EASUR','E    '/
      DATA HD/'SQUAR','ED EU','CLIDI','AN DI','STANC','E    '/ 
      DATA (HS(1,J),J=1,4)/'NEARE','ST NE','IGHBO','UR   '/
      DATA (HS(2,J),J=1,4)/'FURTH','EST N','EIGHB','OUR  '/
      DATA (HS(3,J),J=1,3)/'GROUP',' AVER','AGE  '/
      DATA (HS(4,J),J=1,2)/'MEDIA','N    '/
      DATA (HS(5,J),J=1,2)/'CENTR','OID  '/
      DATA (HS(6,J),J=1,4)/'FLEXI','BLE S','ORTIN','G    '/
      DATA (HS(7,J),J=1,5)/'INCRE','MENTA','L SUM',' OF S','QUARE'/

C                  CLUSTERING METHODS                           
C         'Nearest Neighbour (Single Linkage)          =1              
C         'Furthest Neighbour (Maximum Linkage)        =2              
C         'Group Average (Average Linkage)             =3              
C         'Median                                      =4              
C         'Centriod                                    =5              
C         'Flexible Sorting                            =6              
C         'Incremental Sum of Squares (WARD''S Method) =7             

      FAC = 1.0                                                         
      NM=MAX0(NIND,NATT)                                                    
      N1=1                                                              
      N2=N1+NIND*NATT                                                       
      N3=N2+NM                                                          
      LIM=N3                                                            
      IF(LIM.GT.LIMZ) GOTO 40                                          
      DO 10 I=1,N3                                                      
10      Z(I)=0.0                                                       
15    CALL RCLUH(Z(N1),Z(N2),H1,H2,H3,HH1,HH2,FORM1,NIND,NATT,          
     &           ISU,FAC,X) 
c      write (19,*) 'FORM 1'
c      write (19,*) (FORM1(III),III=1,14)
      ND=(NIND-1)*NIND/2                                                    
      N3=N2+ND                                                          
      N4=N3+NIND                                                          
      N5=N4+NIND                                                          
      N6=N5+2*NIND                                                        
      LIM=N6                                                            
      IF(LIM.GT.LIMZ) GOTO 40                                          
      DO 20 I=N2,N3                                                     
20      Z(I)=0.0                                                        
      CALL DIST(Z(N1),Z(N2),NIND,NATT,FORM3,HD,HH4,X)                 
cccccccccccccc
c      CALL WSIM(Z(N2),Z(N3),NIND)
ccccccccccccccc
      DO 30 I=N3,N6                                                     
30      Z(I)=0.0                                                        
                                                                        
      CALL AGHICL(Z(N2),Z(N3),Z(N4),Z(N5),NIND,BETA,IS)                   
      CALL ALLOC(Z(N5),NIND,NATT,NG,ISU,IS,BETA,H2,HS,X,IDT)
      RETURN
                                                                        
40    WRITE(*,50) LIMZ,LIM                                             
50    FORMAT(1X,'THE DIMENSIONS OF THE Z ARRAY ARE TOO SMALL'/          
     &       1X,'THEY SHOULD BE INCREASED FROM',I6,'TO ',I6)                   
      IFAULT=9
      RETURN                                                              
      END                                                               
                                                                        
                                                                        
      SUBROUTINE RCLUH(RC,TEMP,H1,H2,H3,HH1,HH2,FORM1,NIND,NATT,        
     &                 ISU,FAC,X)                                               
C     Reads GXE matrix, stores as either GXE or EXG, standardizes &/0r     
C     col. corrects if required.                                           
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      DIMENSION RC(*),X(MNIND,MNATT),TEMP(*)                           
      CHARACTER*5 H1(2,3),H2(2,3),H3(2,3),HH1(*),HH2(*),FORM1(*)
                                                                        
      DO 110 I=4,6                                                        
        I1=I-3                                                            
        FORM1(I)=HH1(I1)                                                 
110   CONTINUE
      DO 120 I=10,11                                                      
        I1=I-9                                                            
        FORM1(I)=HH2(I1)                                                 
120   CONTINUE
      DO 130 IR=1,NIND                                                    
        DO 130 IC=1,NATT
          M1=(IR-1)*NATT+IC                                                   
          RC(M1)=X(IR,IC)*FAC                                             
130   CONTINUE
      DO 140 I=1,3                                                       
140     FORM1(I)=H1(1,I)                                                
      DO 150 I=1,3                                                     
        I1=I+11                                                           
        FORM1(I1)=H3(1,I)                                               
150   CONTINUE
      IF(ISU.EQ.1) THEN                                            
        DO 170 I=1,3                                                     
          I1=I+6                                                            
          FORM1(I1)=H2(1,I)                                               
170     CONTINUE
      ELSE
        DO 180 I=1,3                                                       
          I1=I+6                                                            
          FORM1(I1)=H2(2,I)                                               
180     CONTINUE
        CALL STAND(RC,NIND,NATT)                                              
      ENDIF

      DO 190 II=1,NIND
        DO 190 JJ=1,NATT
          TEMP((II-1)*NATT+JJ)=X(II,JJ) 
190   CONTINUE
      RETURN                                                          
      END                                                               

       SUBROUTINE WSIM(D,T,NR)
C
C  WRITES DISSIMILARITY MATRIX.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION D(*),T(*)
C
       WRITE(21,1000)
1000      FORMAT(1H1//1X,20HDISSIMILARITY MATRIX/)
C
	 I1F=NR-1
          DO 10 I1=1,I1F
	     I2S=I1+1
	   I3=0
          DO 11 I2=I2S,NR
	      M1=(I2-1)*(I2-2)/2+I1
	     I3=I3+1
11      T(I3)=D(M1)
           I4S=1
	 I4F=NR-I1
10      WRITE(21,1001)(T(I4),I4=I4S,I4F)
1001      FORMAT(10F10.3)
         RETURN
         END





      SUBROUTINE STAND(A,NIND,NATT)                                         
C     Col. standardizes the A matrix by standard deviation.                
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)                                                    
                                                                        
      DO 220 IC=1,NATT
        N=0.0                                                             
        X=0.0                                                             
        XX=0.0                                                            
        DO 200 IR=1,NIND                                                     
          M=(IR-1)*NATT+IC                                                    
          IF(A(M).LE.-90.0) GOTO 200 
          XX=XX+A(M)**2                                                     
          X=X+A(M)                                                          
          N=N+1.0                                                           
200      CONTINUE                                                        
        XX=SQRT((XX-X**2/N)/(N-1.))                                       
        DO 210 IR=1,NIND                                                     
          M=(IR-1)*NATT+IC                                                    
          IF(A(M).LE.-90.0) GOTO 210                                         
          A(M)=A(M)/XX                                                      
210      CONTINUE                                                        
220    CONTINUE                                                        
      RETURN                                                            
      END                                                               
                                                                        

      SUBROUTINE DEL(A,NIND,NATT)                                           
C     Corrects each row  by mean of row                                    
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)                                                    
      REAL N                                                            
                                                                        
      DO 320 IR=1,NIND                                                     
        X=0.0                                                             
        N=0.0                                                             
        DO 300 IC=1,NATT                                                     
          M=(IR-1)*NATT+IC                                                    
          IF(A(M).LE.-90.0) GOTO 300                                         
          X=X+A(M)                                                          
          N=N+1.0                                                           
300     CONTINUE                                                        
        X=X/N                                                             
        DO 310 IC=1,NATT                                                     
          M=(IR-1)*NATT+IC                                                    
          IF(A(M).LE.-90.0) GOTO 310                                         
          A(M)=A(M)-X                                                       
310     CONTINUE                                                        
320   CONTINUE                                                        
      RETURN                                                            
      END                                                               
 

      SUBROUTINE DIST(RC,DIST1,NIND,NATT,FORM3,HD,HH4,X)              
C     Calculates the dissimilarity matrix as required.                     
C     only squared Euclidean distance implemented;                
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RC(*),DIST1(*)                                          
      CHARACTER*5              FORM3(*),HD(*),HH4(*)                    
                                                                        
      DO 400 I=1,5                                                       
400     FORM3(I)=HH4(I)                                                 
      DO 410 I=1,6                                                       
        I1=I+5                                                            
        FORM3(I1)=HD(I)                                                 
410   CONTINUE
      CALL SED(RC,DIST1,NIND,NATT)                                          
      RETURN                                                            
      END                                                               
 

      SUBROUTINE SED(RC,DIST1,NIND,NATT)                                    
C     Requires an NIND by NATT matrix of data & calculates the SED             
C     dissimilarity matrix between the nind individuals & stores them       
C     in the array DIST.                                                   
C     missing data is denoted by a value .LE.-90.0.                        
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RC(*),DIST1(*)                                          
                                                                        
      I1F=NIND-1                                                          
      DO 510 I1=1,I1F                                                    
        I2S=I1+1                                                          
        DO 510 I2=I2S,NIND                                                   
          M1=(I2-1)*(I2-2)/2+I1                                             
          DIV=NATT                                                            
          DO 500 J=1,NATT                                                      
            M2=(I1-1)*NATT+J                                                    
            M3=(I2-1)*NATT+J                       
	    IF ((RC(M2).LE.-90.0).OR.(RC(M3).LE.-90.0)) THEN    
              DIV=DIV-1.0                                                     
            ELSE
             DIST1(M1)=DIST1(M1)+(RC(M2)-RC(M3))**2                            
            ENDIF 
500        CONTINUE                                                        
          DIST1(M1)=DIST1(M1)/DIV                                         
510    CONTINUE
      RETURN                                                            
      END                                                               
 

      SUBROUTINE AGHICL(D,N,M,MEMB,NT,BET,IS)  
C     Requires the dissimilarity matrix (D), the number of individuals         
C     to be classified (NT), and the nominated BETA value (BET) for        
C     flexible sorting if required.                                        
C      N; Stores the NO. of members in the groups.                         
C      M; Stores group names.                                              
C     Returns memb which stores the hierarchy for further use.             
C     clusters the NT taxa by any one of the nominated agglomerative       
C     hierarchical clustering strategies (IS) & writes the hierarchy        
C     & dissimilarity measure on fusion.                                   
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION N(*),M(*),MEMB(*)
      DIMENSION D(*)               
                                                                        
      DO 600 I=1,NT                                                       
        N(I)=1                                                            
        M(I)=I                                                           
600   CONTINUE
                                                                        
      IF=NT-1                                                           
      DO 610 I=1,IF                                                      
        CALL MIND(D,N,IMIN,JMIN,DMIN,NT)                                  
        M1=(I-1)*2+1                                                      
        M2=M1+1                                                           
        MEMB(M1)=M(IMIN)                                                  
        MEMB(M2)=M(JMIN)                                                  
c	write (19,*) IS
c	write (19,*) MEMB(M1),MEMB(M2)
        NTI=NT+I                                                          
        CALL DISTO(D,N,IMIN,JMIN,NT,BET,IS)                               
        N(IMIN)=N(IMIN)+N(JMIN)                                           
        N(JMIN)=0                                                         
        M(IMIN)=NT+I                                                      
610    CONTINUE                                                        
      RETURN                                                            
      END                                                               
 

      SUBROUTINE MIND(D,N,IMIN,JMIN,DMIN,NT)                            
C     Requires the GP.-GP. dissimilarity matrix D & the array of           
C     number of members in groups (N).                                     
C     returns the names of the two GPS. (IMIN & JMIN) with the             
C     smallest dissimilarity between them & this dissimilarity             
C     (DMIN).                                                              
C     NT = NO. of taxa being clustered                                     
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION N(*)
      DIMENSION D(*)                                             
                                                                        
      DMIN=1.0E+20                                                      
      I1F=NT-1                                                          
      DO 700 I1=1,I1F                                                    
        I2S=I1+1                                                          
        DO 700 I2=I2S,NT                                                   
          IF(N(I1).EQ.0) GOTO 700
          IF(N(I2).EQ.0) GOTO 700 
          M1=(I2-1)*(I2-2)/2+I1                                             
          IF(DMIN.LT.D(M1)) GOTO 700
          DMIN=D(M1)                                                        
          IMIN=I1                                                           
          JMIN=I2                                                           
700   CONTINUE                                                        
      RETURN                                                            
      END                                                               


      SUBROUTINE DISTO(D,N,IMIN,JMIN,NT,BET,IS)                         
C     Calculates the new GP.-GP. distances between the new GP. formed      
C     when GPS. IMIN & JMIN FUSE & all other GPS.                          
C      BET; Nominated value for flexible sorting if required.              
C      IS ; Code for strategy of clustering.                               
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION N(*)
      DIMENSION D(*)                                              
                                                                        
      DO 800 K=1,NT                                                      
        IF ((N(K).NE.0).AND.(K.NE.IMIN).AND.(K.NE.JMIN)) THEN     
          K1=K                                                              
          CALL DISTQ(D,N,IMIN,JMIN,DISS,K1,BET,IS)                          
          IF(IMIN.GT.K) THEN                                             
            M1=(IMIN-1)*(IMIN-2)/2+K                                        
          ELSE
            M1=(K-2)*(K-1)/2+IMIN                                             
          ENDIF
          D(M1)=DISS                                                      
        ENDIF
800   CONTINUE                                                        
      RETURN                                                            
      END                                                               
 
      SUBROUTINE DISTQ(D,N,IMIN,JMIN,DISS,K,BET,IS)                     
C     Calculates new GP.-GP. Dissimilarity between new GP. formed          
C     when GPS. IMIN & JMIN Fuse & other GPS.                              
C       IS = 1; NEAREST NEIGHBOUR (SINGLE LINKAGE)                         
C          = 2; FURTHEST NEIGHBOUR (MAXIMUM LINKAGE)                       
C          = 3; GROUP AVERAGE (AVERAGE LINKAGE)                            
C          = 4; MEDIAN                                                     
C          = 5; CENTROID                                                   
C          = 6; FLEXIBLE SORTING (BET NEEDS TO BE SPECIFIED)               
C          = 7; INCREMENTAL SUM OF SQUARES (BURR'S METHOD OR WARD'S        
C               METHOD)                                                    
C                                                                       
C     See WISHART 1971; BIOMETRICS 25: PP165-70:  For algorithm on         
C     all strategies except flexible sorting:  whence lance &              
C     WILLIAMS 1967; COMPUTER J.: 9 : PP373-80.                            
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION N(*)
      DIMENSION D(*)                                            
      REAL NI,NJ,NK,NM                                               
                                                                        
      NI=N(IMIN)                                                        
      NJ=N(JMIN)                                                        
      NK=N(K)                                                           
      NR=NI+NJ                                                          
      NM=NR+NK                                                          
                                                                        
      IF(IMIN.GT.K) THEN
        M1=(IMIN-1)*(IMIN-2)/2+K                                        
      ELSE
        M1=(K-1)*(K-2)/2+IMIN                                             
      ENDIF                                                          

      IF(JMIN.GT.K) THEN
        M2=(JMIN-1)*(JMIN-2)/2+K                                        
      ELSE                                     
        M2=(K-2)*(K-1)/2+JMIN                                             
      ENDIF
      M3=(JMIN-1)*(JMIN-2)/2+IMIN                                     
                                                                        
      DIK=D(M1)                                                         
      DJK=D(M2)                                                         
      DIKJK=ABS(DIK-DJK)                                                
      DIJ=D(M3)                                                         
                                                                        
      GO TO (1,2,3,4,5,6,7)IS                                           
                                                                        
C  Nearest neighbour                                                    
1      ALA=0.5                                                          
      ALB=0.5                                                           
      BETA=0.0                                                          
      GAM=-0.5                                                          
      GO TO 10                                                          

C  Furthest neighbour                                                   
2      ALA=0.5                                                          
      ALB=0.5                                                           
      BETA=0.0                                                          
      GAM=0.5                                                           
      GO TO 10                                                          
                                                                        
C  Group average                                                        
3      ALA=NI/NR                                                        
      ALB=NJ/NR                                                         
      BETA=0.0                                                          
      GAM=0.0                                                           
      GO TO 10                                                          
                                                                        
C  Median                                                               
4      ALA=0.5                                                          
      ALB=0.5                                                           
c      BETA=0.0                                                          
c      GAM=-0.25                                                         
       GAM=0.0                                                          
       BETA=-0.25                                                         
      GO TO 10                                                          
                                                                        
C  Centroid                                                             
5      ALA=NI/NR                                                        
      ALB=NJ/NR                                                       
      BETA=-1*ALA*ALB                                                      
      GAM=0.0                                                           
      GO TO 10                                                          
                                                                        
C  Flexible sorting                                                     
6      ALA=(1.0-BET)/2.0                                                
      ALB=ALA                                                           
      BETA=BET                                                          
      GAM=0.0                                                           
      GO TO 10                                                          
                                                                        
C  Incremental sum of squares                                           
7      ALA=(NI+NK)/NM                                                   
      ALB=(NJ+NK)/NM                                                    
      BETA=(-NK)/NM                                                       
      GAM=0.0                                                           
                                                                        
C  Calculate DISS Depending on strategy                                 
10      DISS=ALA*DIK+ALB*DJK+BETA*DIJ+GAM*DIKJK                         
      RETURN                                                            
      END                                                               
                                                                        
      
      SUBROUTINE ALLOC(IHIE,NIND,NATT,NG,ISU,IS,BETA,H2,HS,X,IDT)
C     This subroutine takes the hierarchical form and determines a partition
C     of the data set into NG groups

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       INCLUDE 'EMMIX.max'
      DOUBLE PRECISION IHIE(MNIND)
      DIMENSION NAME(MNIND),IDT(MNIND),
     &          NUMB(MNIND),MEMB(MNIND)              
      CHARACTER*5 H2(2,3),HS(7,6),HEAD(9)                                  

      DO 900 I=1,3                                                       
900     HEAD(I)=H2(ISU,I)                                                    
      DO 910 I=1,6                                                       
910    HEAD(I+3)=HS(IS,I)                                                   
      IF (IS.NE.6) BETA=0.0                                             
      CALL FUSE(NAME,NUMB,MEMB,IHIE,NIND,NIND,NG)
      IG1=0       
      DO 920 IF=1,NG                                                     
        DO 920 IG=1,NUMB(IF)                                               
           IG1=IG1+1                                                         
           IG2=MEMB(IG1)                                                     
           IDT(IG2)=IF                                                         
920   CONTINUE
c      CALL CHECK(NG,NIND,IDT)
      RETURN
      END

      SUBROUTINE FUSE(NAME,NUMB,MEMB,IHIE,NIND,NGPS,NGPF)                 
C     Four single dimensional arrays containing                         
C      NAME  :Group names                                               
C      NUMB  :NO. in each group                                         
C      MEMB  :Membership of each group                                  
C      IHIE  :Hierarchy from agglomerative clustering                  
C                                                                       
C      IG1   :Name of first group to fuse at a level                    
C      IG2   :Name of 2nd.    "    "  "   "  "  "                       
C      NGL   :NO. of GPS. before fusion at that level                   
C      NENG  :Name of new GRP.                                          
C      NB1   :Number of memb's. in GPS. before GP. IG1                  
C      NB2   :  "    "    "     "   "     "     "  IG2                  
C      IP1   :POSITION OF GP1 IN NAME & NUMB                            
C      IP2   :  "      "  GP2  "  "   "  "                              
C      NS1   :NO. OF MEMB'S. OF GP1 OF FUSION                           
C      NS2   : "  "   "      "  GP2 "    "                              
C      NNG   : "  "   "      "  NEW GP. AFTER FUSION                    

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION IHIE(*)
      DIMENSION NAME(*),NUMB(*),
     &          MEMB(*)               
                                                                        
      DO 1000 I=1,NIND
        NAME(I)=I                                                         
        NUMB(I)=1                                                         
        MEMB(I)=I                                                         
1000  CONTINUE

      DO 1020 I=NGPS,NGPF+1,-1                                            
        NGP=I                                                             
        NF=NIND-NGP                                                         
        NAGP=NIND+NF                                                        
        NAGPN=NAGP+1                                                      
        IP=NF*2+1                                                         
        IG1=IHIE(IP)                                                      
        IG2=IHIE(IP+1)                                                    
        CALL FUSEL(NAME,NUMB,MEMB,IG1,IG2,NGP,NAGPN,NIND)                   
        NGP=NGP-1                                                         
        NF=NF+1                                                           
1020  CONTINUE                                                          
      RETURN                                                            
      END                                                               


      SUBROUTINE FUSEL (NAME,NUMB,MEMB,IG1,IG2,NGL,NEWG,NIND)             
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NAME(*),NUMB(*),MEMB(*)                           
C     Finding place of grps. to fuse in name & numb array     
      DO 1100 I=1,NGL                                                     
        IF(NAME(I).EQ.IG1)IP1=I                                           
        IF(NAME(I).EQ.IG2)IP2=I                                           
1100  CONTINUE                                                          
      IPL=MAX0(IP1,IP2)                                                 
      IPS=MIN0(IP1,IP2)                                                 
C     Finding place in membership array of members                      
C     of groups to fuse                                                 
      NB1=0                                                             
      DO 1110 I=1,IPL-1                                                   
1110    NB1=NB1+NUMB(I)                                                   
      NB2=0                                                             
      DO 1120 I=1,IPS-1                                                   
1120    NB2=NB2+NUMB(I)                                                   
C     No. of members of gps. to fuse                                    
      IPS=MIN0(IP1,IP2)                                                 
      NS1=NUMB(IPL)                                                     
      NS2=NUMB(IPS)                                                     
      NNG=NS1+NS2                                                       
C     Shifting positions of names & numbers in groups                   
C     to positions after fusion                                         
      DO 1130 I=IPL,NGL-1                                                 
        NAME(I)=NAME(I+1)                                                 
        NUMB(I)=NUMB(I+1)                                                 
1130  CONTINUE                                                 
      DO 1140 I=IPS,NGL-1                                                 
        NAME(I)=NAME(I+1)                                                 
        NUMB(I)=NUMB(I+1)                                                 
1140  CONTINUE
      NAME(NGL-1)=NEWG                                                  
      NUMB(NGL-1)=NNG                                                   
C     Shifting MEMBS. of GPS. to appropriate position of                
C     Memb array                                                        
      CALL SHIFT(MEMB,NB1,NS1,NIND)                                       
      CALL SHIFT(MEMB,NB2,NS2,NIND)                                       
      RETURN                                                            
      END                                                               
 

      SUBROUTINE SHIFT(MEMB,NB,NS,NG)                                   
C     Shifts membership of MEMB according to fusion                     

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION MEMB(*)                                               
                                                                        
      DO 1210 I=1,NS                                                      
        IT=MEMB(NB+1)                                                     
        DO 1200 J=NB+2,NG                                                   
1200       MEMB(J-1)=MEMB(J)                                                 
        MEMB(NG)=IT                                                       
1210  CONTINUE
      RETURN                                                            
      END                                                               
