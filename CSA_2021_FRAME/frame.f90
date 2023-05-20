
! March 03, 2021
! COMPUTATIONAL STRUCTURAL ANALYSIS
! Input file: input.txt
! 2D Frame [m, kN]

! Module is used to define all global variables in the program

MODULE GLOBALVAR
! for GEOMET Subroutine
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)     											 ! All lines starting from A-H and O-Z are in double precision
	INTEGER:: NNODE, NELE, NSEC, NMAT	
	INTEGER, DIMENSION(:), ALLOCATABLE, SAVE:: ISEC, IMAT      						     ! Incidence of Section and Material vectors
	INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE:: IN				 					     ! Incidence of element matrix 
	REAL(2), DIMENSION(:,:), ALLOCATABLE, SAVE:: COORD, CSEC, CMAT					   	 ! (2) double precision
	
! March 10, 2021
! for SCODE Subroutine
	INTEGER:: NDOF																		 ! Total number of degree of freedom
	INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE:: IDOF

! March 15, 2021
! for LOADS Subroutine
	REAL(2), DIMENSION(:,:), ALLOCATABLE, SAVE:: ELOADS									 ! Matrix of element loads
	REAL(2), DIMENSION(:), ALLOCATABLE, SAVE:: VLOADS  									 ! Vector of nodal loads
	REAL(2), DIMENSION(:,:), ALLOCATABLE, SAVE:: PRES									 ! Prestressing
	
! March 23, 2021
! for ASSEMB - MKK Subroutine
	REAL(2), DIMENSION(:,:), ALLOCATABLE, SAVE:: VK											! GLOBAL STIFFNESS MATRIX
	REAL(2), DIMENSION(6,6):: ST															! LOCAL STIFFNESS MATRIX
	REAL(2), DIMENSION(6):: EQFG															! EQUIVALENT NODAL FORCE IN GLOBAL REFERENCE SYSTEM

! March 31, 2021
! for SOLVE Subroutine
	REAL(2), DIMENSION(:), ALLOCATABLE, SAVE:: VDISP										! final vector of nodal loads

END MODULE GLOBALVAR
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! MAIN PROGRAM

PROGRAM FRAME
	USE GLOBALVAR
    ! Local Variable
    CHARACTER(LEN=80):: TITLE										! the variable through which we read the global title of the input
	! Opening and Closing files
    OPEN(11, FILE='input.txt', STATUS='UNKNOWN')
    OPEN(12, FILE='output.txt', STATUS='UNKNOWN')
    OPEN(9, FILE='MKK.tmp', STATUS='UNKNOWN')						! temporary file
    
    READ(11,*)TITLE													! * free format
    
    ! CALL OF SUBROUTINE
    CALL GEOMET
    CALL SCODE      		 ! UPDATE ON MARCH 10, 2021
    CALL LOADS       		 ! UPDATE ON MARCH 15, 2021
    CALL ASSEMB			 	 ! UPDATE ON MARCH 23, 2021
    CALL JOINTS				 ! UPDATE ON APRIL 12, 2021
	CALL SOLVE				 ! UPDATE ON MARCH 31, 2021    
	CALL STRESS				 ! UPDATE ON APRIL 07, 2021
	CALL PLOT				 ! UPDATE ON APRIL 14, 2021
	
    CLOSE(11)
    CLOSE(12)
    CLOSE(9)
    STOP
  
END PROGRAM FRAME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! SUBROUTINES 

SUBROUTINE GEOMET
		USE GLOBALVAR
		INTEGER:: IDNODE, IDSEC, IDMAT, IDELE, I, J, ICODE							! I, J are counter and ICODE is to identify specific error
		CHARACTER(LEN=80)::REC          											! TO READ A LINE FROM THE INPUT.TXT, reading the header of the section        
          ! Read and Write the data in specific format
          	WRITE(12,'(/,"GEOMETRY SECTIONS/MATERIALS")')
          	READ(11,*)REC
          	READ(11,*)NNODE
          		! Verify the input
          		IF (NNODE < 2) THEN
             		WRITE(12,'(/,1X,"DATA ERROR-NODES")')
             		STOP
             		RETURN
          		END IF
          
          	WRITE(12,'(/, 1X, "NO. OF NODES = ", 1X, I3)')NNODE
          	WRITE(12, 101) 'NODE', 'COORD_X', 'COORD_Y'						! Format 101 should be provided after RETURN and before end of subroutine
          
          	ALLOCATE(COORD(NNODE,2))										! Matrix of the coordinate
          
          	DO I=1, NNODE
             	READ (11, *) IDNODE, COORD(IDNODE,:)						! read IDNODE and all columns of the coordinate
             	WRITE (12, 102) IDNODE, COORD(IDNODE,:)
            END DO
             
       		! SECTIONS
       		READ(11,*) REC													! read the next header: SECTIONS
       		READ(11,*) NSEC
       			 ! verify the input
                 IF (NSEC<1) THEN
                 	  WRITE(12,'(/,"DATA ERROR-NSEC")')
                      STOP
                      RETURN
                 END IF
                  
            ALLOCATE(CSEC(NSEC, 4))											! Matrix of the characterstics of the cross-section
                        
       		WRITE(12,'(/, 1X, "NO. OF SECTIONS = ", 1X, I3)')NSEC
       		WRITE(12, 103) 'SEC', 'AREA', 'INERTIA', 'DEPTH', 'SF'
                  
       		DO I = 1, NSEC
          		READ(11,*) IDSEC, CSEC(IDSEC,:)
          		WRITE(12,104) IDSEC, CSEC(IDSEC,:)
          	
          		ICODE = 0			! Verify the quality of the data

       			DO J = 1,3
          			IF (CSEC(IDSEC,J)<=0) ICODE=J 							! Input value at index 'J' is invalid.
       			END DO
       			
       			IF (CSEC(IDSEC,4)<0) ICODE=4								! Shear correction factor should be non-negative
       
       					IF (ICODE > 0 ) THEN
          					WRITE(12,'(/, 1X, "DATA ERROR-SEC.PROP_CSEC(", I2, ",", I2, ")")') IDSEC, ICODE
          					STOP
          					RETURN
          				END IF
          	END DO
          
          	! MATERIAL
          	READ(11,*)REC
          	READ(11,*)NMAT
          		 ! Verify the input
          		 IF (NMAT<1) THEN
                 	  WRITE(12,'(/, 1X, "DATA ERROR-MAT")')
                      STOP
                      RETURN
                 END IF
                 
          	ALLOCATE(CMAT(NMAT,4))											! Martix of the characterstics of the material
          	
          	WRITE(12,'(/, 1X, "NO. OF MATERIALS = ", 1X, I3)')NMAT
       		WRITE(12, 105) 'MAT', 'E', 'v', 'THERMAL', 'SELF_WEIGHT'
          	
          	DO I=1, NMAT
             	READ(11,*) IDMAT, CMAT(IDMAT,:)
             	WRITE(12,106) IDMAT, CMAT(IDMAT,:)
             	
             	ICODE = 0
             	
             	IF(CMAT(IDMAT,1)<=0) ICODE=1											! Young modulus should be positive
             	IF((CMAT(IDMAT,2)<0).OR.(CMAT(IDMAT,2)>0.5)) ICODE=2					! Poisson ratio out of range
             	IF(CMAT(IDMAT,3)<0) ICODE=3												! Thermal coefficient should be positive
             	IF(CMAT(IDMAT,4)<0) ICODE=4												! Self-weight should be positive
             	
             			IF(ICODE>0) THEN
                   			WRITE(12,'(/, 1X, "DATA ERROR-MAT.PROP_CMAT(", I2, ",", I2, ")")') IDMAT, ICODE
                   			STOP
                   			RETURN
             			END IF
          	END DO

			! CONNECTIVITY
			READ(11,*)REC
			READ(11,*)NELE
				! Verify the input
        	 	IF(NELE<1) THEN
                  	  WRITE(12,'(/, 1X, "DATA ERROR-ELEMENT")')
                      STOP
                      RETURN
         		END IF
         
         ALLOCATE(IN(NELE,2))									! Incidence of element matrix
         ALLOCATE(ISEC(NELE))									! Incidence vector of the cross-section
         ALLOCATE(IMAT(NELE))									! Incidence vector of the material
         
         WRITE(12,'(/,1X,"NO. OF ELEM = ", 1X, I3)')NELE
         WRITE(12,107)'ELEM','NODE_1','NODE_2','SECT_YPE','MAT_TYPE'

         DO I=1, NELE
             READ(11,*)IDELE,(IN(IDELE,J),J=1,2),ISEC(IDELE),IMAT(IDELE)							! Alternative way of writing inner do-loop   
             WRITE(12,108)IDELE,(IN(IDELE,J),J=1,2),ISEC(IDELE),IMAT(IDELE)							
         END DO
         
         RETURN													! Return back to the main program
         
        	 ! DEFINITION OF THE FORMAT LIST
		 101 FORMAT (1X, A4, 10X, A12, 5X, A12)
		 102 FORMAT (1X, I4, 10X, F12.2, 5X, F12.2)
		 103 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12, 5X, A12)
		 104 FORMAT (1X, I4, 10X, E12.4, 5X, E12.4, 5X, F12.3, 5X, F12.2)							! Section
		 105 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12, 5X, A12)
		 106 FORMAT (1X, I4, 10X, E12.4, 5X, F12.2, 5X, E12.4, 5X, F12.2)							! Material 
		 107 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12, 5X, A12)
		 108 FORMAT (1X, I4, 10X, I12, 5X, I12, 5X, I12, 5X, I12)

END SUBROUTINE GEOMET

!!! March 10, 2021

SUBROUTINE SCODE
		USE GLOBALVAR
        INTEGER::I, J, NRES, NLINK, IDNODE, JDIR, NMAST, NSLAVE
        CHARACTER(LEN=80)::REC     											 				! READ THE HEADER OF DIFFERENT SECTION  
           ! INITIALIZATION
           ALLOCATE(IDOF(NNODE,3))											 				! MATRIX OF DEGREE OF FREEDOM
           IDOF(:,:)=0
           WRITE(12,'(/,"RESTRAINTS/LINKS")')
           ! READ/WRITE RESTRAINTS
           READ(11,*)REC
           READ(11,*)NRES           										 				! NUMBER OF RESTRAINTS
           WRITE(12,'(/, 1X, "NO. OF RESTRAINTS = ", 1X, I3)')NRES
           WRITE(12, 101) 'NODE', 'DIRECTION'
           
           DO I=1,NRES
              READ(11,*)IDNODE,JDIR
              WRITE(12,102)IDNODE,JDIR
              IDOF(IDNODE,JDIR)=-1
           END DO
           
           ! READING/WRITING LINKS
           READ(11,*)REC
           READ(11,*)NLINK
           			! Verify the input
           			IF(NLINK<0)THEN
              			 WRITE(12,'("DATA ERROR-NO. OF LINK")')
               			 STOP
               			 RETURN
          		    END IF
           
           WRITE(12,'(/,1X,"NO. OF LINKS = ", 1X, I3)')NLINK
           
           IF (NLINK>0) THEN																! It is possible to have no links
              WRITE(12,103)'MASTER','SLAVE','DIR'
              DO I=1, NLINK
                 READ(11,*) NMAST, NSLAVE, JDIR
                 WRITE(12,104) NMAST, NSLAVE, JDIR
                 IDOF(NSLAVE,JDIR)=NMAST
              END DO
           END IF
           
           ! Numbering of DOFs
		   NDOF=0
		   WRITE(12,'(/,1X,"DOF MATRIX")')
		   WRITE(12,105) 'Node', 'Ux', 'Uy', 'Rotz'
		   
		   EXT: DO I = 1, NNODE
     			INT: DO J = 1, 3
            			IF(IDOF(I,J)<0) CYCLE INT											! CYCLE jumps to the next iteration if restrained (-1)
            			IF(IDOF(I,J)==0) THEN
                			NDOF=NDOF+1
                			IDOF(I,J)=NDOF
            			ELSE
                			IDNODE=IDOF(I,J)												! IDOF(I,J) is the number of the master node
                			IDOF(I,J)=IDOF(IDNODE,J)
            			END IF
          			 END DO INT
              		 WRITE(12,106) I, IDOF(I,:)
          		END DO EXT
          		RETURN
          		
           ! DEFINITION OF THE FORMAT
		 101 FORMAT (1X, A4, 10X, A12)
		 102 FORMAT (1X, I4, 10X, I12)
		 103 FORMAT (1X, A4, 10X, A12, 5X, A12)
		 104 FORMAT (1X, I4, 10X, I12, 5X, I12)
		 105 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12)
		 106 FORMAT (1X, I4, 10X, I12, 5X, I12, 5X, I12)
		 
END SUBROUTINE SCODE

! MARCH 15, 2021

SUBROUTINE LOADS
		USE GLOBALVAR
		! Local variable declaration
        INTEGER:: I, NNL, IDNODE, JDIR, JDOF, NEL, IDELE
        INTEGER:: NEP, NTEND, IDTEND					 	 				! For prestressing load
        REAL(2):: PFORCE, ELOAD(5)                        					! ELOAD is NOT ELOADS
        REAL(2):: PRE(4), FP, E1, EM, E2									! PRES is NOT PRE
        CHARACTER(LEN = 80):: REC
           
           WRITE(12,'(/,"NODAL AND ELEMENT LOADS")')
           
           ALLOCATE(VLOADS(NDOF))											! Vector of Nodal loads in global reference system
           ALLOCATE(ELOADS(NELE, 5))										! Matrix of Element loads: [Px, Py1, Py2, DTop, DBot]
           ALLOCATE(PRES(NELE,4))											! Matrix of Prestressing loads
           
           VLOADS(:)=0.              										! Initialization (DOUBLE PRECISION)
           ELOADS(:,:)=0.
           PRES(:,:)=0.
           
           READ(11,*)REC
           READ(11,*)NNL             										! Total number of Nodal loads
           			! Verify the input
           			IF (NNL<0) THEN
              			WRITE(12,'("DATA ERROR-NO. OF LOAD")')
              			STOP
              			RETURN
           			END IF
           
           WRITE(12,'(/,1X,"NO. OF LOADS = ", 1X, I3)')NNL
           
           IF (NNL>0) THEN
              WRITE(12,101) 'NODE', 'DIR', 'FORCE'
              DO I = 1, NNL
                   READ (11,*) IDNODE, JDIR, PFORCE
                   IF((IDNODE<1).OR.(IDNODE>NNODE))THEN
                         WRITE(12,'("NODAL LOADS APPLIED ON INCORRECT POSITION.")')
                         STOP
                         RETURN
                   END IF
                   
                   IF((JDIR<1).OR.(JDIR>3))THEN
                         WRITE(12,'("NODAL LOADS APPLIED ALONG INCORRECT DIRECTION")')
                         STOP
                         RETURN
                   END IF
                   
                   WRITE(12, 102) IDNODE, JDIR, PFORCE
                   JDOF=IDOF(IDNODE, JDIR)
                   VLOADS(JDOF)=VLOADS(JDOF) + PFORCE          				! It is possible to have two nodal loads at a point
              END DO
              
              WRITE(12,'(/, 1X, "VECTOR OF NODAL LOADS")')
              WRITE(12,103) 'DOF', 'FORCE'
              DO I = 1, NDOF
                   WRITE(12, 104) I, VLOADS(I)
              END DO
           END IF
           
          ! ELEMENT LOADS THAT ARE DEFINED IN A LOCAL REFERENCE SYSTEM
          
           READ(11,*)REC
           READ(11,*)NEL													! Total number of element loads
           			! Verify the input
	                 IF (NEL<0) THEN
                  		  WRITE(12,'("DATA ERROR-NO. ELEMENT LOAD")')
                    	  STOP
                          RETURN
                     END IF
                 
                 WRITE(12,'(/, 1X, "NO. OF ELEMENT LOADS = ", 1X, I3)')NEL
                 IF(NEL>0)THEN
                       WRITE(12, 105)'ELEM', 'Px', 'Py1', 'Py2', 'DTtop', 'DTbot'
                       DO I = 1, NEL
                             READ(11,*) IDELE, ELOAD    		         						! ALL CHARACTERSTIC OF ELOADS ARE SAVED
                             IF ((IDELE<1).OR.(IDELE>NELE)) THEN
                                   WRITE(12, '("ID NUMBER OF LOADED ELEMENT IS INCORRECT")')
                                   STOP
                                   RETURN
                             END IF
                        	 WRITE(12, 106) IDELE, ELOAD										! ELOAD(5)
                        	 ELOADS(IDELE,1:5) = ELOADS(IDELE,1:5) + ELOAD(1:5)					! ELOADS(5)
                        END DO
                 END IF     
                 
            ! UPDATE ON MARCH 17, 2021
       		! PRESTRESSING
       		READ(11,*)REC
       		READ(11,*)NEP													! Number of prestressed element
       				 ! Verify the input
	                 IF (NEP<0) THEN
                  		  WRITE(12,'("ERROR-NO. OF PRESTRESSED ELEMENT")')
                    	  STOP
                          RETURN
                     END IF
       		
       			WRITE(12,'(/,1X,"NO. OF PRESTRESSED ELEMENT = ", 1X, I3)')NEP
       			IF (NEP>0) THEN
       				DO I=1,NEP
       					READ(11,*)IDELE, NTEND
       					IF ((IDELE<1).OR.(IDELE>NELE)) THEN
       						WRITE(12, '("ID NUMBER OF PRESTRESSED ELEMENT IS INCORRECT")')
       						STOP
       						RETURN
       					END IF
       				
       				FP = 0.
       				E1 = 0.
       				EM = 0.
       				E2 = 0.
       					
       				WRITE(12, '(1X, "ELEMENT = ", 1X, I3)')IDELE
       				WRITE(12, 107)'TENDON', 'P', 'e1', 'em', 'e2'
       				
       				DO J = 1, NTEND
       					READ(11,*) IDTEND, PRE(1:4)
       					WRITE(12, 108) IDTEND, PRE(1:4)
       					FP = FP + PRE(1)
       					E1 = E1 + PRE(1)*PRE(2)
       					EM = EM + PRE(1)*PRE(3)
       					E2 = E2 + PRE(1)*PRE(4)
       				END DO
       				
       				PRES(IDELE,1)=FP
       				PRES(IDELE,2)=E1/FP
       				PRES(IDELE,3)=EM/FP       				
       				PRES(IDELE,4)=E2/FP
       		
       				WRITE(12,109)'EQUIVALENT', FP, E1/FP, EM/FP, E2/FP
       				END DO
       			END IF
       		 RETURN
        
	   ! DEFINITION OF THE FORMAT LIST
		 101 FORMAT (1X, A4, 10X, A12, 5X, A12)
		 102 FORMAT (1X, I4, 10X, I12, 5X, F12.2)
		 103 FORMAT (1X, A4, 10X, A12)
		 104 FORMAT (1X, I4, 10X, F12.2)
		 105 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12, 5X, A12, 5X, A12)
		 106 FORMAT (1X, I4, 10X, F12.2, 5X, F12.2, 5X, F12.2, 5X, F12.2, 5X, F12.2)  
		 107 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12, 5X, A12)  
		 108 FORMAT (1X, I4, 10X, F12.2, 5X, F12.2, 5X, F12.2, 5X, F12.2)   
		 109 FORMAT (1X, A4, 10X, F12.2, 5X, F12.2, 5X, F12.2, 5X, F12.2) 
       
END SUBROUTINE LOADS																! NOW EXPECT THE COMPLETE DESCRIPTION OF THE FRAME STRUCTURE
       
! UPDATE ON MARCH 23, 2021

SUBROUTINE ASSEMB
		USE GLOBALVAR
		INTEGER:: I, J, NE, NCODE(6)
		
		ALLOCATE(VK(NDOF, NDOF))													! GLOBAL STIFFNESS MATRIX
		! INITIALIZATION
		VK(:,:)=0.
		DO NE = 1, NELE																! NO. OF ELEMENT
			CALL MKK(NE)															! THE OUTPUT OF MKK ARE EQFG(6) AND ST(6,6) LOCAL STIFFNESS MATRIX IN GLOBAL SYSTEM
			N1 = IN(NE,1)															! IN(NELE,2): it contains the nodal points of an element
			N2 = IN(NE,2)															! Identify the row element
			NCODE(1:3)=IDOF(N1,1:3)
			NCODE(4:6)=IDOF(N2,1:3)
			
			EXT: DO I=1,6
				IC = NCODE(I)
				IF(IC<0) CYCLE EXT			! CYCLE if constrained since the global stiffness matrix should not be modified
				VLOADS(IC)= VLOADS(IC) + EQFG(I) 			! CHECK EQFG(I)
				INT: DO J=1,6
					JC = NCODE(J)
					IF(JC<0) CYCLE INT		! CYCLE if constrained
					VK(IC,JC) = VK(IC,JC) + ST(I,J)
				END DO INT
			END DO EXT
		END DO
		
		! PRINT OF VK AND VLOADS
		WRITE(12,'(/)')
		WRITE(12,101)'DOF', 'FORCE'
		WRITE(12,102)(I, VLOADS(I), I=1,NDOF)
		WRITE(12,'(/,1X,"STIFFNESS MATRIX")')					! Shorter alternative to write DO loop
			
			DO J = 1,NDOF
				WRITE(12,'(/,1X,"ROW",T18,"COLUMN",I3)')J		! the Output of the stiffness matrix is displayed as a table
				DO I = 1,NDOF
					WRITE(12,103) I, VK(I,J)
				END DO
			END DO
			
		RETURN
		
		! DEFINITION OF THE FORMAT LIST
		 101 FORMAT (1X, A4, 10X, A12)
		 102 FORMAT (1X, I4, 10X, E12.4)
		 103 FORMAT (1X, I4, 10X, E12.4)
		
END SUBROUTINE ASSEMB

! MKK

SUBROUTINE MKK(NE)
		USE GLOBALVAR
		! Local variable declaration
		INTEGER:: I, J, N1, N2, ISC, IMT
		REAL(2):: DX, DY, AL, CA, SA, AA, AJ, AH, CHI, EMOD, POISSON, ALPHA, WEIGHT, GMOD, PHI, PX, PY1, PY2, DTBOT, DTTOP, & 
					& FP, E1, EM, E2, TETA1, TETA2, CURV, WW, TMP(6,6), T(6,6), EQF(6)			! EQF (Local) is NOT EQFG (Global)
		
			! Geometric properties
			N1 = IN(NE, 1)
			N2 = IN(NE, 2)
			DX = COORD(N2,1) - COORD(N1,1)
			DY = COORD(N2,2) - COORD(N1,2)
			AL = SQRT(DX**2 + DY**2)			! LENGTH OF THE ELEMENT (NE)
			CA = DX/AL							! COSINE
			SA = DY/AL							! SINE
			
			! Sections and Materials     
			ISC = ISEC(NE)						! INCIDENCE OF CROSS-SECTION
			IMT = IMAT(NE)						! INCIDENCE OF MATERIAL
			AA = CSEC(ISC,1)					! AREA
			AJ = CSEC(ISC,2)					! INERTIA
			AH = CSEC(ISC,3)					! DEPTH
			CHI = CSEC(ISC,4)					! SHEAR FACTOR
			EMOD = CMAT(IMT,1)					! YOUNG MODULES
			POISSON = CMAT(IMT,2)				! POISSON RATIO
			ALPHA = CMAT(IMT,3)					! THERMAL COEFFICIENT
			WEIGHT = CMAT(IMT,4)				! SELF-WEIGHT
			GMOD = EMOD/(2.*(1 + POISSON))				! SHEAR MODULUS
			PHI = 12.*CHI*EMOD*AJ/(AL**2.*GMOD*AA)		! SHEAR-DEFORMABILITY COEFFICIENT .(DOUBLE PRECISION)
			
			! Local Stiffness Matrix (tiangular lower matrix)
				ST(:,:) = 0.							! local Stiffness Matrix initialization
				ST(1,1) = EMOD*(AA/AL)
				ST(2,2) = 12.*EMOD*AJ/((1+PHI)*AL**3.)
				ST(3,2) = 6.*EMOD*AJ/((1+PHI)*AL**2.)
				ST(3,3) = (4.+PHI)*EMOD*AJ/((1.+PHI)*AL)
				ST(4,1) = -ST(1,1)
				ST(4,4) = ST(1,1)
				ST(5,2) = -ST(2,2)
				ST(5,3) = -ST(3,2)
				ST(5,5) = ST(2,2)
				ST(6,2) = ST(3,2)
				ST(6,3) = (2.-PHI)*EMOD*AJ/((1.+PHI)*AL)
				ST(6,5) = -ST(3,2)
				ST(6,6) = ST(3,3)
			! Fill in the other upper triangular matrix exploiting symmetry
				DO I=1,6
					DO J=I+1,6							! the diagonal elements are excluded
						ST(I,J) = ST(J,I)
					END DO
				END DO
				
				WRITE(9,*) ST(:,:)						! Write the local stiffness matrix on temporary file
				
			! Transformation matrix
				T(:,:) = 0.
				T(1,1) = CA
				T(1,2) = -SA
				T(2,1) = SA
				T(2,2) = CA
				T(3,3) = 1.
				T(4,4) = CA
				T(4,5) = -SA
				T(5,4) = SA
				T(5,5) = CA
				T(6,6) = 1.
				
				WRITE(9,*) T(:,:)
			! Global Stiffness matrix, K = TKT^T 
				TMP = MATMUL(T,ST)						! MATMUL: Embedded FORTRAN function to multiply two matrices
				ST = MATMUL(TMP, TRANSPOSE(T))			! TRANSPOSE: Embedded FORTRAN function
				
			! Equivalent nodal forces
				EQF(:) = 0.
			! Distributed loads
				PX = ELOADS(NE,1)
				PY1 = ELOADS(NE,2)
				PY2 = ELOADS(NE,3)
				DTTOP = ELOADS(NE,4)
				DTBOT = ELOADS(NE,5)
			! Selfweight influence
				PG = WEIGHT*AA							! Weight of the element
				PX = PX - PG*SA							! Contribution of the self-weight is added to each load 
				PY1 = PY1 - PG*CA
				PY2 = PY2 - PG*CA
				QQ = PY2 - PY1
			! Prestressing forces
				FP = PRES(NE,1)
				E1 = PRES(NE,2)
				EM = PRES(NE,3)
				E2 = PRES(NE,4)
			! Contribution of prestressing
				TETA1 = -1./AL*(3.*E1 + E2 - 4.*EM)
				TETA2 = 1./AL*(E1 + 3.*E2 - 4.*EM)
				CURV = 4.*(E1 + E2 - 2.*EM)/(AL**2)
				WW = FP*CURV
				PY1 = PY1 + WW
				PY2 = PY2 + WW
			! Longitudinal loads (in x-direction)
				EQF(1) = EQF(1) + PX*AL/2.
				EQF(4) = EQF(4) + PX*AL/2.
			! Transveral loads (CONSTANT)
				EQF(2) = EQF(2) + PY1*AL/2.
				EQF(3) = EQF(3) + PY1*AL**2./12.
				EQF(5) = EQF(5) + PY1*AL/2.
				EQF(6) = EQF(6) - PY1*AL**2./12.
			! Transversal loads (LINEAR)
				COEF_A = AL/(12.*EMOD*AJ)
				COEF_B = CHI/(AL*AA*GMOD)
				COEF_C = COEF_B/(COEF_A + COEF_B)
				EQF(2) = EQF(2) + 9./60.*QQ*AL*(1.+COEF_C/9.)
				EQF(3) = EQF(3) + 1./30.*QQ*AL**2.*(1.+COEF_C/4.)
				EQF(5) = EQF(5) + 21./60.*QQ*AL*(1.-COEF_C/21.)
				EQF(6) = EQF(6) - 1./20.*QQ*AL**2.*(1.-COEF_C/6.)
			! Thermal loads (elongation)
				DELTA_TM = (DTTOP + DTBOT)/2.
				FH = EMOD*AA*ALPHA*DELTA_TM
				EQF(1) = EQF(1) - FH 
				EQF(4) = EQF(4) + FH
			! Transveral loads (flexural behaviour)
				DELTA_TD = (DTTOP - DTBOT)/2.
				FM = 2.*ALPHA*DELTA_TD*EMOD*AJ/AH
				EQF(3) = EQF(3) + FM
				EQF(6) = EQF(6) - FM
				
				WRITE(9,*) EQF(:)
			! The subsequent modification does NOT affect internal stress
				EQF(1) = EQF(1) + FP
				EQF(2) = EQF(2) + FP*TETA1
				EQF(3) = EQF(3) - FP*E1
				EQF(4) = EQF(4) - FP
				EQF(5) = EQF(5) - FP*TETA2
				EQF(6) = EQF(6) + FP*E2
			! Global equivalent nodal load, f=Tf'
				EQFG = MATMUL(T, EQF)
		RETURN
		
END SUBROUTINE MKK

! Subroutine SOLVE

SUBROUTINE JOINTS
		USE GLOBALVAR
		INTEGER:: I, J, JDIR
		INTEGER:: NELR, NELL, NTDIS, NRDIS
		INTEGER:: IDNODE, N1, N2
		REAL(2):: STIFF, ADIS, VKJJ
		CHARACTER(LEN=80):: REC
		WRITE(12,'(/,T1,"ELASTIC RESTRAINTS-LINKS")')
		
		! ELASTIC RESTRAINTS
		READ(11,*)REC
		READ(11,*)NELR
		IF(NELR<0)THEN
			WRITE(12,'("DATA ERROR-NO. OF ELASTIC RESTRAINTS")')
		END IF
		
		WRITE(12,'(1X,"NO. OF ELASTIC RESTRAINTS = ", 1X, I3)')NELR
		IF(NELR>0) THEN
			WRITE(12,101)'NODE','DIRECTION','STIFFNESS'
			DO I=1,NELR
				READ(11,*)IDNODE, JDIR, STIFF
				WRITE(12,102)IDNODE, JDIR, STIFF
				JDOF = IDOF(IDNODE,JDIR)
				VK(JDIR,JDOF) = VK(JDOF, JDOF) + STIFF
			END DO
		END IF
		
		! ELASTIC LINKS
		READ(11,*)REC
		READ(11,*)NELL
		IF(NELL<0)THEN
			WRITE(12,'("DATA ERROR-ELASTIC LINKS")')
		END IF
		
		WRITE(12,'(/,1X,"NO. OF ELASTIC LINKS = ", 1X, I3)')NELL
		IF (NELL>0)THEN
			WRITE(12,103) 'NODE1', 'NODE2', 'DIRECTION', 'STIFFNESS'
			DO I=1, NELL
				READ(11,*) N1,N2,JDIR,STIFF
				WRITE(12,104) N1,N2,JDIR,STIFF
				JDOF1 = IDOF(N1, JDIR)
				JDOF2 = IDOF(N2, JDIR)
				VK(JDOF1, JDOF1) = VK(JDOF1, JDOF1) + STIFF
				VK(JDOF2, JDOF2) = VK(JDOF2, JDOF2) + STIFF
				VK(JDOF1, JDOF2) = VK(JDOF1, JDOF2) - STIFF
				VK(JDOF2, JDOF1) = VK(JDOF1, JDOF2)
			END DO
		END IF

		! IMPOSED DISP(ABSOLUTE)
		READ(11,*)REC
		READ(11,*)NTDIS
		IF(NTDIS<0)THEN
			WRITE(12,'("DATA ERROR-NO. OF IMPOSED DISPLACEMENT")')
		END IF
		
		WRITE(12,'(/,1X,"NO. OF IMPOSED DISPLACEMENT = ", 1X, I3)')NTDIS
		IF(NTDIS>0) THEN 
			WRITE(12, 105) 'NODE', 'DIRECTION', 'DISPALCEMENT'
			DO I=1, NTDIS
				READ(11,*) IDNODE, JDIR, ADIS
				WRITE(12, 106) IDNODE, JDIR, ADIS
				JDOF = IDOF(IDNODE, JDIR)
				VKJJ = VK(JDOF, JDOF)
				VK(JDOF,1:NDOF) = 0.
				VLOADS(1:NDOF) = VLOADS(1:NDOF) - VK(1:NDOF, JDOF)*ADIS							! CHECK NDOF Vs DOF
				VK(1:NDOF, JDOF) = 0.
				VK(JDOF, JDOF) = VKJJ
				VLOADS(JDOF) = VKJJ*ADIS
			END DO
		END IF
		
		! IMPOSED DISPLACEMENT(RELATIVE)
		READ(11,*)REC
		READ(11,*)NRDIS
		IF (NRDIS<0)THEN
			WRITE(12, '("DATA ERROR-NO. OF RELATIVE IMPOSED DISPLACEMENT")')
		END IF
		
		WRITE(12,'(/1X,"NO. OF RELATIVE DISPLACEMENT = ", 1X, I3)')NRDIS
		IF (NRDIS>0)THEN
			WRITE(12,107)'NODE1', 'NODE2', 'DIRECTION', 'DISPLACEMENT'
			DO I = 1,NRDIS
				READ(11,*) N1, N2, JDIR, ADIS
				WRITE(12, 108) N1, N2, JDIR, ADIS
				IF ((N1<1).OR.(N1>NNODE))THEN
					WRITE(12,'("DATA ERROR-INCORRECT NODES")')
				END IF
				IF((N2<1).OR.(N2>NNODE))THEN
					WRITE(12,'("DATA ERROR-INCORRECT NODES")')
				END IF
				IF (N1==N2)THEN
					WRITE(12,'("DATA ERROR-INCORRECT NODES")')
				END IF
				IF ((JDIR<1).OR.(JDIR>3))THEN
					WRITE(12,'("DATA ERROR-INCORRECT DIRECTION")')
				END IF
				JDOF1 = IDOF(N1, JDIR)
				JDOF2 = IDOF(N2, JDIR)
				IF ((JDOF1==-1).OR.(JDOF2==-1))THEN
					WRITE(12,'("DATA ERROR-")')
				END IF
				
				VKJJ = VK(JDOF1, JDOF1) + VK(JDOF2, JDOF2)
				DO K = 1, NDOF
					VK(JDOF1,K) = VK(JDOF1,K) + VK(JDOF2,K)
				END DO
				
				VLOADS(JDOF1) = VLOADS(JDOF1) + VLOADS(JDOF2)
				DO K = 1, NDOF														! CHECK I Vs J
					VK(JDOF2, K) = 0
					VK(K, JDOF1) = VK(K,JDOF1) + VK(K,JDOF2)
					VLOADS(K) = VLOADS(K) - VK(K,JDOF2)*ADIS
					VK(K, JDOF2) = 0.
				END DO
				! APPLY ALL MODIFICATION TO THE STIFFENSS MATRIX AND VECTOR OF LOADS
				VK(JDOF1, JDOF1) = VK(JDOF1, JDOF1) + VKJJ
				VK(JDOF1, JDOF2) = -VKJJ
				VK(JDOF2, JDOF1) = -VKJJ
				VK(JDOF2, JDOF2) = VKJJ
				VLOADS(JDOF1) = VLOADS(JDOF1) - VKJJ*ADIS
				VLOADS(JDOF2) = VKJJ*ADIS
			END DO
		END IF 
		
		! INCLINED CONSTRAINTS
		
		RETURN
		! FORMAT LISTS
		101 FORMAT (1X, A4, 10X, A12, 5X, A12)
		102 FORMAT (1X, I4, 10X, I12, 5X, F12.5)
		103 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12)
		104 FORMAT (1X, I4, 10X, I12, 5X, I12, 5X, F12.5)
		105 FORMAT (1X, A4, 10X, A12, 5X, A12)
		106 FORMAT (1X, I4, 10X, I12, 5X, F12.2)  
		107 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12)  
		108 FORMAT (1X, I4, 10X, I12, 5X, I12, 5X, F12.5)
		
END SUBROUTINE JOINTS

! Sobroutine Solve

SUBROUTINE SOLVE
		USE GLOBALVAR
		INTEGER:: I, J, K, JDOF
		REAL(2):: DNOD(3), DMX(3), DMN(3), PVAL				! We are also looking for the maximum and minimum displacement
		
		ALLOCATE(VDISP(NDOF))
		! Initialization of the vector of the constant terms
		! Load vector - > Displacement vector			
		VDISP = VLOADS					! VDISP is intially the modified vector of known terms
	    ! Gaussian elimination method
			DO I=1,NDOF					! I is the counter of the dominant 
				DO K=I+1,NDOF
					VK(I,K) = VK(I,K)/VK(I,I)
				END DO
				VDISP(I) = VDISP(I)/VK(I,I)
				VK(I,I) = 1.
				DO J=I+1,NDOF
					DO K = I+1,NDOF
						VK(J,K) = VK(J,K)-VK(J,I)*VK(I,K)
					END DO
					VDISP(J) = VDISP(J)-VK(J,I)*VDISP(I)
					VK(J,I) = 0.								! The lower part of the matrix is zero
				END DO
			END DO
			
		! Backward substitution
			DO I=NDOF,1,-1
				DO K=I+1,NDOF
					VDISP(I) = VDISP(I)-VK(I,K)*VDISP(K)
				END DO
			END DO												! The vector of nodal load is defined by now
		! Write the output
		WRITE(12,'(/,1X, "SUMMARY OF RESULTS")')
		WRITE(12,'(1X, "NODAL DISPLACEMENTS")')
		WRITE(12, 101) 'NODE', 'Ux', 'Uy', 'Rotz'
		
		DNOD = 0.
		DMX = -1.E6												! Even a dispplacement of zero is considered maximum
		DMN = 1.E6												! Even the node is restrained
		
		DO I=1, NNODE
			DO J=1, 3
				IF (IDOF(I,J)<0)THEN							! The node is restrained thus the displacement is zero
					DNOD(J)=0.
				ELSE
					JDOF=IDOF(I,J)
					DNOD(J)=VDISP(JDOF)
				END IF
				PVAL = DMX(J)
				DMX(J) = MAX(PVAL, DNOD(J))
				PVAL = DMN(J)
				DMN(J) = MIN(PVAL, DNOD(J))
			END DO
			WRITE(12, 102) I, (DNOD(J), J=1,3)
		END DO
					
		WRITE(12, 103) 'MAX', (DMX(J), J=1,3)
		WRITE(12, 104) 'MIN', (DMN(J), J=1,3)
		WRITE(12,'(/,1X,"DISPLACEMENTS AT DEGREE OF FREEDOM")')
		WRITE(12, 105) 'DOF', 'DISP'
		
		DO I=1, NDOF
			WRITE(12, 106)I, VDISP(I)
		END DO
		
		RETURN
		
		! DEFINATION OF FORMAT LISTS
		101 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12)
		102 FORMAT (1X, I4, 10X, E12.5, 5X, E12.5, 5X, E12.5)
		103 FORMAT (1X, A4, 10X, E12.5, 5X, E12.5, 5X, E12.5)
		104 FORMAT (1X, A4, 10X, E12.5, 5X, E12.5, 5X, E12.5)
		105 FORMAT (1X, A4, 10X, A12)
		106 FORMAT (1X, I4, 10X, E12.5)  
		
END SUBROUTINE SOLVE

! Subroutine STRESS

SUBROUTINE STRESS
		USE GLOBALVAR
		INTEGER:: I, J, JDOF, NCODE(6)
		REAL(2):: QS(6), T(6,6), EQF(6), QF(6), QSL(6)
		
		REWIND(9)												! REwIND THE TEMPORARY FILE TO EXTRACT LOCAL STIFFNESS, TRANSFORMATION MATRIX AND EQF
		
		WRITE(12,'(/,1X,"FORCES AT NODES")')
		WRITE(12,101) 'ELEM', 'N1', 'T1', 'M1', 'N2', 'T2', 'M2'
		
		DO NE=1,NELE
			N1 = IN(NE,1)
			N2 = IN(NE,2)
			NCODE(1:3) = IDOF(N1,1:3)
			NCODE(4:6) = IDOF(N2,1:3)
			READ(9,*)ST
			READ(9,*)T
			READ(9,*)EQF
			
			DO I=1,6
				JDOF=NCODE(I)
				IF(JDOF==-1)THEN
					QS(I)=0.
				ELSE
					QS(I)=VDISP(JDOF)
				END IF
			END DO
			
			QSL = MATMUL(TRANSPOSE(T),QS)
			
			DO I=1,6
				QF(I)=-EQF(I)
				DO J=1,6
					QF(I) = QF(I) + ST(I,J)*QSL(J)
				END DO
			END DO
			
			WRITE(12,102) NE, (QF(J), J=1,6)
		END DO
		RETURN
		
		! DEFINATION OF FORMAT LISTS
		101 FORMAT (1X, A4, 10X, A12, 5X, A12, 5X, A12, 5X, A12, 5X, A12, 5X, A12)
		102 FORMAT (1X, I4, 10X, E12.5, 5X, E12.5, 5X, E12.5, 5X, E12.5, 5X, E12.5, 5X, E12.5)
	
END SUBROUTINE STRESS

! SUBROUTINE PLOT

SUBROUTINE PLOT
		USE GLOBALVAR
		IMPLICIT DOUBLE PRECISION (A-H,O-Z)
		INTEGER:: I, J, JDOF, NCODE(6)
		REAL(2):: QS(6), T(6,6), EQF(6), QF(6), QSL(6)
		! Rewind of the temporary file
		REWIND(9)
		! Opening file 33 (input data for PLOT.EXE)
		OPEN (33,FILE='PLOT.DAT',STATUS='UNKNOWN')
		! Writing of coordinates and incidence
		WRITE (33,*)NNODE,NELE
		DO I=1,NNODE
			WRITE (33,*) (COORD(I,J),J=1,2)
		END DO
		DO NE=1,NELE
			WRITE (33,*) (IN(NE,J),J=1,2)
		END DO
		! Loop for all the elements
		DO NE=1,NELE
			N1 = IN(NE,1)
			N2 = IN(NE,2)
			NCODE(1:3) = IDOF(N1,1:3)
			NCODE(4:6) = IDOF(N2,1:3)
			! Reading of the stiffness matrix
			READ (9,*) ST
			! Reading of the transformation matrix
			READ (9,*) T
			! Reading of the equivalent nodal forces
			READ (9,*) EQF
			! Nodal displacements in the global ref. syst.
			DO I=1,6
				JDOF = NCODE(I)
				IF(JDOF==-1) THEN
					QS(I) = 0.
				ELSE
					QS(I) = VDISP(JDOF)
				ENDIF
			END DO
			
			! Nodal displacements in the local ref. syst.
			QSL = MATMUL(TRANSPOSE(T),QS)
			! Forces at end nodes
			DO I=1,6
				QF(I) = -EQF(I)
				DO J=1,6
					QF(I) = QF(I)+ST(I,J)*QSL(J)
				END DO
			END DO
			
			! Displacements and forces along the length of the element
			DX = COORD(N2,1)-COORD(N1,1)
			DY = COORD(N2,2)-COORD(N1,2)
			AL = SQRT(DX**2.+DY**2.)
			CA = DX/AL
			SA = DY/AL
			ISC = ISEC(NE)
			IMT = IMAT(NE)
			AA = CSEC(ISC,1)
			WEIGHT = CMAT (IMT,4)
			PX = ELOADS(NE,1)
			PY1 = ELOADS(NE,2)
			PY2 = ELOADS(NE,3)
			PG = AA*WEIGHT
			PX = PX-PG*SA
			PY1 = PY1-PG*CA
			PY2 = PY2-PG*CA
			FP = PRES(NE,1)
			E1 = PRES(NE,2)
			EM = PRES(NE,3)
			E2 = PRES(NE,4)
			CURV = 4.*(E1+E2-2.*EM)/(AL**2.)
			WW = FP*CURV
			PY1 = PY1+WW
			PY2 = PY2+WW
			
				DO I=0,20
					XX=AL*DBLE(I)/20.
					PY=PY1+(PY2-PY1)*XX/AL
					UU=U(QSL,XX,AL)
					VV=V(QSL,XX,AL)
					SSN=-QF(1)-PX*XX
					SST=QF(2)+(PY1+PY)*XX/2.
					SSM=-QF(3)+QF(2)*XX+PY1*XX**2./2.+(PY-PY1)*XX**2./6.
					WRITE (33,*) UU,VV,SSN,SST,SSM
				END DO
		END DO
		CLOSE (33)
		RETURN
END SUBROUTINE PLOT

!====================================
! HERMITE FUNCTIONS
!====================================
DOUBLE PRECISION FUNCTION F(I,AL,X)
			IMPLICIT DOUBLE PRECISION (A-H,O-Z)
			CSI = X/AL
			SELECT CASE (I)
				CASE(1)
					F = 1.-CSI
				CASE(4)
					F = CSI
				CASE(2)
					F = 1.-3.*CSI**2.+2.*CSI**3.
				CASE(3)
					F = AL*(CSI-2.*CSI**2.+CSI**3.)
				CASE(5)
					F = 3.*CSI**2.-2.*CSI**3.
				CASE(6)
					F = AL*(CSI**3.-CSI**2.)
			END SELECT
END FUNCTION F

DOUBLE PRECISION FUNCTION U(QSL,X,AL)
			IMPLICIT DOUBLE PRECISION (A-H,O-Z)
			DIMENSION :: QSL(6)
			U = 0.
			DO I=1,4,3
				U = U+F(I,AL,X)*QSL(I)
			END DO
		RETURN
END FUNCTION U

DOUBLE PRECISION FUNCTION V(QSL,X,AL)
			IMPLICIT DOUBLE PRECISION (A-H,O-Z)
			DIMENSION :: QSL(6)
			V = 0.
			DO I=2,3
				V = V+F(I,AL,X)*QSL(I)
			END DO
			DO I=5,6
				V = V+F(I,AL,X)*QSL(I)
			END DO
		RETURN
END FUNCTION V
				
			
		
 
