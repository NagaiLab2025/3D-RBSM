!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     �R�����q�a�r�l��̓v���O����
!     3D RBSM PROGRAM - SIMULATION PART
!           VERSION 2.1 (DIRECT METHOD BY PARDISO)  10/2012
!           VERSION  MSM045   FOR BIAXIAL          10/2004
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
!       ��@  �@    �@���@��@�@�G�@��
!    COPYRIGHT          KOHEI NAGAI
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
MODULE PARAM_CONST

    IMPLICIT NONE
    REAL,PARAMETER::PI=3.14159265
    INTEGER,PARAMETER::DP=KIND(1.0D0)
    INTEGER,PARAMETER::DR_K=selected_real_KIND(14)
    
END MODULE PARAM_CONST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INCLUDE 'mkl_pardiso.f90'
PROGRAM MAIN
      
USE PARAM_CONST
USE mkl_pardiso

IMPLICIT NONE

INTEGER MAX,IREC
INTEGER NNODE,NELE,NPHASE,NXFIX,NYFIX,NZFIX,NSFIX,NXDISP,NYDISP,NZDISP,IFSTEP,MTYPE,IBTYPE,ISTEPBI,ITMAT
INTEGER NNN,NNZ
INTEGER,ALLOCATABLE::NNOPHA(:),INOPHA(:,:),NPHELEM(:),IPHELEM(:,:),IBANE(:),IELPHA(:,:),KIND0(:),&
    IXFIX(:),IYFIX(:),IZFIX(:),ISFIX(:),IXDISP(:),IYDISP(:),IZDISP(:),IXDISPI(:),IYDISPI(:),IZDISPI(:)
INTEGER,ALLOCATABLE::ISKY(:),IK11(:),IK12(:)
REAL(KIND=DP),ALLOCATABLE::COD(:,:),XDISP(:),YDISP(:),ZDISP(:),PRO(:,:)
REAL(KIND=DP),ALLOCATABLE::AREA(:),HH(:,:),CPHASE(:,:),XYZ(:,:)
REAL(KIND=DP),ALLOCATABLE::TEN(:),STDEF(:),STDEF1(:),STDEFBI(:),STDEF1BI(:),BANE(:,:)
    
REAL(KIND=DP),ALLOCATABLE::XK11(:,:),XK12(:,:)
INTEGER,ALLOCATABLE::NK11(:,:),NK12(:,:)!,LPV0(:),LPM0(:)
REAL(DR_K),ALLOCATABLE::RDU(:)
    
INTEGER IXX,I1,I2  !  , Location                     
REAL XX !, XXpor

MAX=100000


WRITE(*,*)'==================================================='
WRITE(*,*)'         3D RBSM SIMULATION ANALYSIS               '
WRITE(*,*)'                         VER. 2.1(DIRECT METHOD)   '
WRITE(*,*)'==================================================='

OPEN(100,FILE='INDATA1.INFO',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)	     
IREC=1
READ(100,REC=IREC) NNODE   !NNODE
IREC=IREC+1
READ(100,REC=IREC) NELE    !NELE
IREC=IREC+1
READ(100,REC=IREC) NPHASE  !NPHASE
CLOSE(100)
      
!!!!!!!!READ DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE(IXFIX(MAX),IYFIX(MAX),IZFIX(MAX),ISFIX(MAX),IXDISP(MAX),IYDISP(MAX),IZDISP(MAX),&
     IXDISPI(MAX),IYDISPI(MAX),IZDISPI(MAX))
ALLOCATE(NNOPHA(NPHASE),INOPHA(NPHASE,50),IELPHA(NPHASE,2),IBANE(NPHASE),NPHELEM(NELE),&
     IPHELEM(NELE,70),KIND0(NELE))
ALLOCATE(COD(NNODE,3),XDISP(MAX),YDISP(MAX),ZDISP(MAX),PRO(13,3))
ALLOCATE(RDU(NELE*6))

!     READ FROM INDATA.INFO
OPEN(100,FILE='INDATA1.INFO',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)	     
IREC=1
IREC=IREC+1
IREC=IREC+1
!     �ߓ_���W
!     COORDINATE OF NODE
DO I1=1,NNODE
    DO I2=1,3
        IREC=IREC+1
	    READ(100,REC=IREC) XX 
	    COD(I1,I2)=XX    
  END DO
END DO
!     �ʍ\���ߓ_�ԍ�
!     NODE NUMBER WHICH COMPOSES FACE
DO I1=1,NPHASE
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   NNOPHA(I1)=IXX
   DO I2=1,NNOPHA(I1)
      IREC=IREC+1
      READ(100,REC=IREC) IXX
      INOPHA(I1,I2)=IXX   
   END DO
END DO
!     �v�f�\���ʔԍ�
!     FACE NUMBER WHICH COMPOSES ELEMENT
DO I1=1,NELE
   IREC=IREC+1
   READ(100,REC=IREC)IXX
   NPHELEM(I1)=IXX
   DO I2=1,NPHELEM(I1)
      IREC=IREC+1
      READ(100,REC=IREC)IXX
      IPHELEM(I1,I2)=IXX        
   END DO
END DO
!     �ʍ\���v�f�ԍ�
!     ELEMENT NUMBER WHICH COMPOSES FACE
DO I1=1,NPHASE
   IREC=IREC+1
   READ(100,REC=IREC)IXX
   IELPHA(I1,1)=IXX
   IREC=IREC+1
   READ(100,REC=IREC)IXX
   IELPHA(I1,2)=IXX     
   IF(IELPHA(I1,2).EQ.0)THEN
      IBANE(I1)=1
   ELSE
      IBANE(I1)=0
   END IF
END DO
!     �v�f��ޔԍ�
!     ELEMENT KIND
DO I1=1,NELE
   IREC=IREC+1
   READ(100,REC=IREC)IXX	  
   KIND0(I1)=IXX 
 
END DO
!     �Œ�v�f����
!     FIXED ELEMENT CONDITION
IREC=IREC+1
READ(100,REC=IREC) NXFIX
IREC=IREC+1
READ(100,REC=IREC) NYFIX
IREC=IREC+1
READ(100,REC=IREC) NZFIX
IREC=IREC+1
READ(100,REC=IREC) NSFIX

print *, "NXFIX=", NXFIX
print *, "NYFIX=", NYFIX
print *, "NZFIX=", NZFIX
print *, "NSFIX=", NSFIX
	    

DO I1=1,NXFIX
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   IXFIX(I1)=IXX
END DO
DO I1=1,NYFIX
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   IYFIX(I1)=IXX
END DO
DO I1=1,NZFIX
   IREC=IREC+1
   READ(100,REC=IREC)IXX
   IZFIX(I1)=IXX
END DO
DO I1=1,NSFIX
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   ISFIX(I1)=IXX
END DO
!     FORCE DISPLACEMENT CONDITION
!     IXDISPI(I1)=   1: TOP SIDE
!                    2: BOTTOM SIDE
!                    3: RIGHT SIDE
!                    4: LEFT SIDE
!                    5: FRONT SIDE
!                    6: BACK SIDE
IREC=IREC+1
READ(100,REC=IREC) NXDISP
IREC=IREC+1
READ(100,REC=IREC) NYDISP
IREC=IREC+1
READ(100,REC=IREC) NZDISP

DO I1=1,NXDISP
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   IXDISP(I1)=IXX
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   IXDISPI(I1)=IXX
   IREC=IREC+1
   READ(100,REC=IREC) XX
   XDISP(I1)=XX
END DO
DO I1=1,NYDISP
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   IYDISP(I1)=IXX
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   IYDISPI(I1)=IXX
   IREC=IREC+1
   READ(100,REC=IREC) XX
   YDISP(I1)=XX
END DO
DO I1=1,NZDISP
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   IZDISP(I1)=IXX
   IREC=IREC+1
   READ(100,REC=IREC) IXX
   IZDISPI(I1)=IXX
   IREC=IREC+1
   READ(100,REC=IREC) XX
   ZDISP(I1)=XX   
END DO
!     FINAL STEP
IREC=IREC+1
READ(100,REC=IREC) IFSTEP
IFSTEP=4000 !4000
!!!!!!!!!!!!!!!!!!
!IFSTEP=100 !! test
!!!!!!!!!!!!!!!!!!
!     FINAL STEP OF APPLING DISPLACEMENT ON RIGHT AND LEFT BOUNDARY
IREC=IREC+1
READ(100,REC=IREC) ISTEPBI
ISTEPBI=4000
CLOSE(100)

!!!!!!READ FROM INDATA2.INFO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN(200,FILE='INDATA2.INFO',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)
!     MODEL TYPE
IREC=1	
READ(200,REC=IREC) MTYPE
IREC=IREC+1
READ(200,REC=IREC) IBTYPE
!     MATERIAL PROPERTIES
IREC=IREC+1
READ(200,REC=IREC) XX
PRO(1,1)=XX  !!MORTAL ELLASTIC MODULAS
IREC=IREC+1 
READ(200,REC=IREC) XX
PRO(1,2)=XX  !!MORTAL POISSON'S RATIO
IREC=IREC+1
READ(200,REC=IREC) XX
PRO(1,3)=XX  !!MORTAL TELSILE STRENGTH
	
IF(MTYPE.EQ.21)THEN
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(6,1)=XX       !! REBAR ELLASTIC MODULAS
    IREC=IREC+1  
    READ(200,REC=IREC) XX
    PRO(6,2)=XX       !! REBAR POISSON'S RATIO
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(4,3)=XX       !! INTERFACE TENSILE STRENGTH
ELSEIF(MTYPE.EQ.22)THEN
    IREC=IREC+1
    READ(200,REC=IREC) XX	
    PRO(6,1)=XX       !! REBAR ELLASTIC MODULAS
    IREC=IREC+1  
    READ(200,REC=IREC) XX
    PRO(6,2)=XX       !! REBAR POISSON'S RATIO
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(4,1)=XX     !! INTERFACE  C
    IREC=IREC+1 
    READ(200,REC=IREC) XX
    PRO(4,2)=XX       !! INTERFACE FAI
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(4,3)=XX       !! INTERFACE STRENGTH
ELSEIF(MTYPE.EQ.31.OR.MTYPE.EQ.41)THEN
    IREC=IREC+1
    READ(200,REC=IREC) XX	
    PRO(2,1)=XX       !!  REBAR ELLASTIC MODULAS
    IREC=IREC+1  
    READ(200,REC=IREC) XX
    PRO(2,2)=XX       !!  REBAR POISSON'S RATIO
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,3)=XX       !!  REBAR TENSILE STRENGTH
    IREC=IREC+1 
    READ(200,REC=IREC) XX
    PRO(4,3)=XX	      !!  INTERFACE TENSILE STRENGTH
ELSEIF(MTYPE.EQ.32.OR.MTYPE.EQ.42)THEN
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,1)=XX       !!  REBAR ELLASTIC MODULAS
    IREC=IREC+1  
    READ(200,REC=IREC) XX
    PRO(2,2)=XX       !!  REBAR POISSON'S RATIO
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,3)=XX       !!  REBAR TENSILE STRENGTH
    IREC=IREC+1 
    READ(200,REC=IREC) XX
    PRO(4,1)=XX       !!  INTERFACE C
    IREC=IREC+1	 
    READ(200,REC=IREC) XX
    PRO(4,2)=XX       !!  INTERFACE FAI
    IREC=IREC+1	 
    READ(200,REC=IREC) XX
    PRO(4,3)=XX	      !!  INTERFACE TENSILE STRENGTH
ELSEIF(MTYPE.EQ.61)THEN     !!����CUTT OFF MODEL
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,1)=XX       !!  REBAR ELLASTIC MODULAS
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,2)=XX       !!  REBAR POISSON'S RATIO
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,3)=XX       !!  REBAR TENSILE STRENGTH
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(4,3)=XX       !!  INTERFACE TENSILR STRENGTH
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(6,1)=XX       !!  AGGREGATE ELLASTIC MODULAS
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(6,2)=XX       !!  AGGREGATE POISSON'S RATIO
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(6,3)=XX       !!  AGGREGATE TENSILE STRENGTH
ELSEIF(MTYPE.EQ.62)THEN     !!����C-FAI MODEL
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,1)=XX       !!  REBAR ELLASTIC MODULAS
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,2)=XX       !!  REBAR POISSON'S RATIO
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(2,3)=XX       !!  REBAR TENSILE STRENGTH
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(4,1)=XX       !!  INTERFACE C
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(4,2)=XX       !!  INTERFACE FAI
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(4,3)=XX       !!  INTERFACE TENSILE STRENGTH
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(6,1)=XX       !!  AGGREGATE ELLASTIC MODULAS
    IREC=IREC+1
    READ(200,REC=IREC) XX
    PRO(6,2)=XX       !!  AGGREGATE POISSON'S RATIO
    IREC=IREC+1 
    READ(200,REC=IREC) XX
    PRO(6,3)=XX       !!  AGGREGATE TENSILE STRENGTH
END IF		

CLOSE(200)
!   ANY VALUE IS OK FOR BOUNDARY ELEMENT
PRO(3,1) = 1.d0
PRO(3,2) = 1.d0
PRO(3,3) = 1.d0
!TENSILE STRENGTH OF INTERFACE IS HALF
PRO(7,1) = PRO(6,1)
PRO(7,2) = PRO(6,2)
!PRO(7,3)=0.5*PRO(1,3)

    !change22112017
    print *, "PRO(1,1)=",PRO(1,1)
    print *, "PRO(1,2)=",PRO(1,2)
    print *, "PRO(1,3)=",PRO(1,3)
    print *, "PRO(6,1)=",PRO(6,1)
   print *, "PRO(6,2)=",PRO(6,2)

  


!Reduce modulus of elasticity and tensile strength for the weak zone please check first with the concrete model
!PRO(13,1)=0.5*PRO(1,1) !reduce modulus of elasticity weak zone
!PRO(13,2)=PRO(1,2)
!PRO(13,3)=0.5*PRO(1,3) !reduce tensile strength weak zone    
    
    
ALLOCATE(REAL(KIND=DP)::AREA(NPHASE),HH(NPHASE,3),CPHASE(NPHASE,3),XYZ(NELE,3))
!--------INDATA PART----------------------------------------------------------------------------------
CALL INDATA(NNODE,NELE,NPHASE,COD,NNOPHA,INOPHA,NPHELEM,IPHELEM,IELPHA,XYZ,AREA,HH,IBANE,CPHASE,ITMAT)  
      
!--------ANALYSIS PART--------------------------------------------------------------------------------			
print *, "     SUBROUTINE ANALYSIS"
      
ALLOCATE(TEN(NPHASE),BANE(NPHASE,2))         
!     CALCULATION OF SPRING STIFFNESS
CALL SPRINGSTIF(NPHASE,IBANE,IELPHA,PRO,HH,BANE,KIND0,XYZ)
      
!     GIVE VARIATION OF SPRING STIFFNESS AND STRENGTH
CALL RANDOMSPRING(MTYPE,TEN,NPHASE,IBANE,IELPHA,PRO,BANE,KIND0,XYZ)
      
ALLOCATE(ISKY(6*NELE),STDEF(6*NELE),STDEF1(6*NELE),STDEFBI(6*NELE),STDEF1BI(6*NELE))
ISKY    = 0
STDEF   = 0.D0
STDEF1  = 0.D0
STDEFBI = 0.D0
STDEF1BI= 0.D0
!     MAKE INDEX
CALL SKY(ISKY,STDEF,STDEF1,STDEFBI,STDEF1BI,NXFIX,NYFIX,NZFIX,NSFIX,IXFIX,IYFIX,IZFIX,ISFIX,&
        NXDISP,NYDISP,NZDISP,IXDISP,IYDISP,IZDISP,IXDISPI,IYDISPI,IZDISPI,XDISP,YDISP,ZDISP,NELE)

!ALLOCATE(LPV0(6*NELE),LPM0(6*NELE)) 
!      MAKE LPV AND LPM FOR SMS-V2 USING ISKY
!      CALL MKLPVLPM(ISKY,LPV0,LPM0,NELE)
    
ALLOCATE(IK11(6*NELE),IK12(6*NELE))
ALLOCATE(NK11(6*NELE,ITMAT))
ALLOCATE(NK12(6*NELE,18))
IK11=0
IK12=0
NK11=0
NK12=0
!     MAKE K11 AND K12
CALL KLOCATION(NELE,ITMAT,ISKY,NPHELEM,IPHELEM,IELPHA,IBANE,IK11,IK12,NK11,NK12)
      
ALLOCATE(XK11(6*NELE,ITMAT))
XK11=0.D0
ALLOCATE(XK12(6*NELE,18))
XK12=0.D0
!     MAKE GLOBAL MATRIX
CALL MAKEMATRIX(NPHASE,NELE,COD,AREA,BANE,IBANE,NNOPHA,INOPHA,IELPHA,ISKY,NK11,XK11,&
                NK12,XK12,IK11,IK12,XYZ,HH,CPHASE,ITMAT)

!       CALCULATION OF NNZ AND NNN FOR SOLVER
CALL CALNZND(NNZ,NELE,NK11,XK11,IK11,ITMAT,NNN,ISKY)

!n     NK11,XK11,NK12,XK12�̓��e���t�@�C���ɏ�������
!n     WRITE NK11, XK11, NK12 AND XK12
!       CALL RECNKXK(NK11,XK11,NK12,XK12,IK11,IK12,NELE,ITMAT)

DEALLOCATE(IXFIX,IYFIX,IZFIX,ISFIX)!,IXDISP,IYDISP,IZDISP,IXDISPI,IYDISPI,IZDISPI,&
           ! XDISP,YDISP,ZDISP)

!--------SIMULATION PART--------------------------------------------------------------------------------------
CALL SIMULATION(NNODE,NELE,NPHASE,COD,NNOPHA,INOPHA,NPHELEM,IPHELEM,IELPHA,KIND0,&
                   ! NXFIX,NYFIX,NZFIX,NSFIX,IXFIX,IYFIX,IZFIX,ISFIX,
                    NXDISP,NYDISP,NZDISP,IXDISP,IYDISP,IZDISP,IXDISPI,IYDISPI,IZDISPI,XDISP,YDISP,ZDISP,&
                    IFSTEP,XYZ,AREA,HH,IBANE,CPHASE,MTYPE,IBTYPE,&
                    PRO,ITMAT,BANE,ISKY,STDEF,STDEF1,STDEFBI,STDEF1BI,IK11,IK12,NNN,TEN,NNZ,ISTEPBI,&
                    NK11,NK12,XK11,XK12)

!11  IF(MODE.EQ.2)THEN
!        CALL RELOAD()
!        CALL RESIMULATION()    
!    ENDIF

STOP 
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE INDATA(NNODE,NELE,NPHASE,COD,NNOPHA,INOPHA,NPHELEM,&
       IPHELEM,IELPHA,XYZ,AREA,HH,IBANE,CPHASE,ITMAT)
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::NNODE,NELE,NPHASE
    INTEGER,INTENT(IN)::NNOPHA(:),INOPHA(:,:),NPHELEM(:),IPHELEM(:,:),IELPHA(:,:)
    INTEGER,INTENT(INOUT)::IBANE(:)
    INTEGER,INTENT(OUT)::ITMAT
    REAL(KIND=DP),INTENT(IN)::COD(:,:)
    REAL(KIND=DP),INTENT(OUT)::XYZ(:,:),AREA(:),HH(:,:),CPHASE(:,:)
!     �e�ʂ̏d�S�̌v�Z
!     CALCULATION OF CENTER OF GRAVITY OF FACE
    CALL PHASECENTER(COD,NPHASE,NNOPHA,INOPHA,CPHASE)     
!     �e�v�f�̏d�S�̌v�Z
!     CALCULATION OF CENTER OF GRAVITY OF ELEMENT
    CALL ELEMCENTER(NELE,NPHELEM,IPHELEM,CPHASE,XYZ)
!     �e�ʂ̖ʐς̌v�Z
!     CALCULATION OF AREA OF FACE
    CALL PHASEAREA(NPHASE,NNOPHA,INOPHA,CPHASE,COD,AREA,IBANE,HH)
!     �e�ʂ���d�S�ւ̐����̋����̌v�Z
!     CALCULATION OF LENGTH OF PEPENDICULAR LINE
    CALL PERPENDICULAR(NPHASE,IBANE,CPHASE,INOPHA,NNOPHA,COD,IELPHA,XYZ,HH)
!     K11N,K11P�̔z��̑傫�������肷��ITMAT�̌v�Z
!     CALCULATION OF ITMAT WHICH DECIDE THE SIZE OF K11N AND K11P
    CALL TMAT(NELE,ITMAT,NPHELEM)

    RETURN
  END SUBROUTINE INDATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PHASECENTER(COD,NPHASE,NNOPHA,INOPHA,CPHASE)
    IMPLICIT NONE

    REAL(KIND=DP),INTENT(IN)::COD(:,:)
    INTEGER,INTENT(IN)::NNOPHA(:),INOPHA(:,:),NPHASE
    REAL(KIND=DP),INTENT(OUT)::CPHASE(:,:)
      
    INTEGER I1,I2,IA,IA1
    REAL(KIND=DP) CX,CY,CZ

!     �e�ʂ̒��S���W�����߂�
!     CALCULATION OF CENTER OF GRAVITY OF FACE
    DO I1=1,NPHASE
       CX=0.0
       CY=0.0
       CZ=0.0
       IA1=NNOPHA(I1)
       DO I2=1,IA1
          IA=INOPHA(I1,I2)
          CX=CX+COD(IA,1)
          CY=CY+COD(IA,2)
          CZ=CZ+COD(IA,3)
       END DO
       CX=CX/IA1
       CY=CY/IA1
       CZ=CZ/IA1
       CPHASE(I1,1)=CX
       CPHASE(I1,2)=CY
       CPHASE(I1,3)=CZ
    END DO

    RETURN
  END SUBROUTINE PHASECENTER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ELEMCENTER(NELE,NPHELEM,IPHELEM,CPHASE,XYZ)
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::NPHELEM(:),IPHELEM(:,:),NELE
    REAL(KIND=DP),INTENT(IN)::CPHASE(:,:)
    REAL(KIND=DP),INTENT(OUT)::XYZ(:,:)
      
    INTEGER I1,I2,IA,IB
    REAL(KIND=DP) XX,YY,ZZ

    DO I1=1,NELE
       XX=0
       YY=0
       ZZ=0
       IA=NPHELEM(I1)
       DO I2=1,IA
          IB=IPHELEM(I1,I2)
          XX=XX+CPHASE(IB,1)
          YY=YY+CPHASE(IB,2)
          ZZ=ZZ+CPHASE(IB,3)
       END DO
       XX=XX/IA
       YY=YY/IA
       ZZ=ZZ/IA
       
       XYZ(I1,1)=XX
       XYZ(I1,2)=YY
       XYZ(I1,3)=ZZ
    END DO

    RETURN
  END SUBROUTINE ELEMCENTER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PHASEAREA(NPHASE,NNOPHA,INOPHA,CPHASE,COD,AREA,IBANE,HH)    
    IMPLICIT NONE
    
    INTEGER,INTENT(IN)::NNOPHA(:),INOPHA(:,:),NPHASE
    REAL(KIND=DP),INTENT(IN)::CPHASE(:,:),COD(:,:)
    REAL(KIND=DP),INTENT(OUT)::HH(:,:),AREA(:)
    INTEGER,INTENT(INOUT)::IBANE(:)
      
    INTEGER I1,I2,IA,IB1,IB2
    REAL(KIND=DP) CX,CY,CZ,SS,X1,X2,Y1,Y2,Z1,Z2,A1,A2,A3,B1,B2,B3,DS

!     �e�ʂ̖ʐς��v�Z
!     CALCULATION OF AREA OF FACE
    DO I1=1,NPHASE
       CX=CPHASE(I1,1)
       CY=CPHASE(I1,2)
       CZ=CPHASE(I1,3)
       SS=0.0
       IA=NNOPHA(I1)
!     �e�ӂƒ��S���W�����񂾎O�p�`�̖ʐς����߁C��������
!     SUM UP TRIANGLE AREAS
       DO I2=1,IA
          IF(I2.NE.IA)THEN
             IB1=INOPHA(I1,I2)
             X1=COD(IB1,1)
             Y1=COD(IB1,2)
             Z1=COD(IB1,3)
             IB2=INOPHA(I1,I2+1)
             X2=COD(IB2,1)
             Y2=COD(IB2,2)
             Z2=COD(IB2,3)
          ELSE
             IB1=INOPHA(I1,I2)
             X1=COD(IB1,1)
             Y1=COD(IB1,2)
             Z1=COD(IB1,3)
             IB2=INOPHA(I1,1)
             X2=COD(IB2,1)
             Y2=COD(IB2,2)
             Z2=COD(IB2,3)
          END IF
!       �x�N�g��
          A1=CX-X1
          A2=CY-Y1
          A3=CZ-Z1
          B1=CX-X2
          B2=CY-Y2
          B3=CZ-Z2
!       �O�p�`�̌v�Z
!       AREA OF TRIANGLE
          DS=SQRT((A2*B3-A3*B2)**2+(A3*B1-A1*B3)**2+(A1*B2-A2*B1)**2)/2
          SS=SS+DS
       END DO
       AREA(I1)=SS

!     ���͌v�Z�p�ɁC�ʂ̂��~�Ɖ��肵���ꍇ�̒��a�����߂�
!     CALCULATION OF DIAMETER WHEN THE FACE SHAPE IS ASSUMED AS CIRCLE
       HH(I1,3)=SQRT(4*SS/3.14159265)
       IF(SS.EQ.0)THEN
	  IBANE(I1)=1
       END IF
    END DO

    RETURN
  END SUBROUTINE PHASEAREA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PERPENDICULAR(NPHASE,IBANE,CPHASE,INOPHA,NNOPHA,COD,IELPHA,XYZ,HH)    
    IMPLICIT NONE

    INTEGER,INTENT(IN)::IBANE(:),NNOPHA(:),INOPHA(:,:),IELPHA(:,:),NPHASE
    REAL(KIND=DP),INTENT(IN)::COD(:,:),XYZ(:,:),CPHASE(:,:)
    REAL(KIND=DP),INTENT(OUT)::HH(:,:)
      
    INTEGER I1,I2,IA,IB1,IB2
    REAL(KIND=DP) A,B,C,D,CX,CY,CZ,X1,X2,Y1,Y2,Z1,Z2,XN,YN,ZN,X0,Y0,Z0,HH2,HH1
!   �e�ʂ���d�S�ւ̐����̒��������߂�
!   CALCULATION OF LENGTH OF PERPENDICULAR FROM THE FACE TO CENTER OF GRAVITY OF ELEMENT
    DO 10 I1=1,NPHASE
!      �o�l��ݒ肵�Ȃ��ʂɑ΂��Ă͌v�Z���Ȃ�
!      NON-SPRING FACE IS NOT CALCULATED
       IF(IBANE(I1).EQ.0)THEN
!       �ʂ̕����������߂�
!       EQUATION OF FACE
          CX=CPHASE(I1,1)   
          CY=CPHASE(I1,2)  
          CZ=CPHASE(I1,3)        
          IA=NNOPHA(I1)
!         �e�ʂƒ��S���W�����񂾎O�p�`�𗘗p���Ė@���x�N�g�������߂�
!         CALCULATION OF NORMAL VECTOR
          DO I2=1,IA
             IF(I2.NE.IA)THEN
                IB1=INOPHA(I1,I2)
                X1=COD(IB1,1)
                Y1=COD(IB1,2)
                Z1=COD(IB1,3)
                IB2=INOPHA(I1,I2+1)
                X2=COD(IB2,1)
                Y2=COD(IB2,2)
                Z2=COD(IB2,3)
             ELSE
                IB1=INOPHA(I1,I2)
                X1=COD(IB1,1)
                Y1=COD(IB1,2)
                Z1=COD(IB1,3)
                IB2=INOPHA(I1,1)
                X2=COD(IB2,1)
                Y2=COD(IB2,2)
                Z2=COD(IB2,3)
             END IF
!            �x�N�g��
!            VECTOR
             X1=CX-X1
             Y1=CY-Y1
             Z1=CZ-Z1
             X2=CX-X2
             Y2=CY-Y2
             Z2=CZ-Z2
!            �@���x�N�g��
!            NORMAL VECTOR
             XN=Y1*Z2-Z1*Y2
             YN=Z1*X2-X1*Z2
             ZN=X1*Y2-Y1*X2
!            AX + BY + CZ +D =0�@�Ƃ����
!            IF AX + BY + CZ +D =0
             A=XN
             B=YN
             C=ZN
             D=-A*CX-B*CY-C*CZ
!            �d�S����̋��������߂�
!            CALCULATION OF LENGTH FROM CENTER OF GRAVITY OF ELEMENT
             HH2=SQRT(A**2+B**2+C**2)
             IF(HH2.NE.0)THEN
                IA=IELPHA(I1,1)
                X0=XYZ(IA,1)
                Y0=XYZ(IA,2)
                Z0=XYZ(IA,3)
                HH1=A*X0+B*Y0+C*Z0+D
                IF(HH1.LT.0)THEN
                   HH1=-HH1
                END IF
                HH(I1,1)=HH1/HH2
                IA=IELPHA(I1,2)
                X0=XYZ(IA,1)
                Y0=XYZ(IA,2)
                Z0=XYZ(IA,3)   
                HH1=A*X0+B*Y0+C*Z0+D
                IF(HH1.LT.0)THEN
                   HH1=-HH1
                END IF
                HH(I1,2)=HH1/HH2
                GOTO 10
             END IF
          END DO
       END IF
10  END DO

    RETURN
  END SUBROUTINE PERPENDICULAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE TMAT(NELE,ITMAT,NPHELEM)
    IMPLICIT NONE
      
    INTEGER,INTENT(IN)::NELE,NPHELEM(:)
    INTEGER,INTENT(OUT)::ITMAT
      
    INTEGER I1,JJ

!   ITMAT�̌v�Z
!   CALCULATION OF ITMAT
    ITMAT=0
!$omp parallel
!$omp do
    DO I1=1,NELE
       JJ=NPHELEM(I1)
       IF(ITMAT.LT.JJ)THEN
          ITMAT=JJ
       END IF
    END DO
!$omp end do
!$omp enda parallel

!     ITMAT���s��̑傫���ɂ���
!     CONVERT TO MATRIX SIZE
    ITMAT=(ITMAT+1)*6

    RETURN
  END SUBROUTINE TMAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SPRINGSTIF(NPHASE,IBANE,IELPHA,PRO,HH,BANE,KIND0,XYZ)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::IBANE(:),IELPHA(:,:),KIND0(:),NPHASE
    REAL(KIND=DP),INTENT(INOUT)::PRO(:,:)
    REAL(KIND=DP),INTENT(IN)::HH(:,:),XYZ(:,:)
    REAL(KIND=DP),INTENT(OUT)::BANE(:,:) 
      
    INTEGER I1,M1,M2,K1,K2
    REAL(KIND=DP) X1,X2,X3,E,P,H,XX,YY,ZZ

!     �|�A�\����ƒe���W���̕␳������
!     MODIFY POISSON RATIO AND ELASTIC MODULUS
!     �e���W���͍Ō��1.05�{����
!     FINALLY, ELASTIC MODULUS IS MULTIPLIED BY 1.05
      DO 20 I1=1,2 

      IF(I1.NE.2)THEN
	  X1=PRO(I1,2)
      X2=-24.8*X1*X1*X1*X1+31.9*X1*X1*X1-16.4*X1*X1+4.28*X1
	  X3=-33.7*X2*X2*X2*X2+17.0*X2*X2*X2-4.13*X2*X2+0.327*X2+1.0
        IF(X3.LE.0.05)THEN
		X3=0.05
		END IF
	  PRO(I1,2)=X2
	  PRO(I1,1)=X3*PRO(I1,1)     
      ELSE
      PRO(2,2)=0
      END IF

20    CONTINUE
    I1=5
    X1=PRO(I1,2)
    X2=-24.8*X1*X1*X1*X1+31.9*X1*X1*X1-16.4*X1*X1+4.28*X1
    X3=-33.7*X2*X2*X2*X2+17.0*X2*X2*X2-4.13*X2*X2+0.327*X2+1.0
        IF(X3.LE.0.05)THEN
        X3=0.05
        END IF
    PRO(I1,2)=X2
    PRO(I1,1)=X3*PRO(I1,1)  
    
      I1=6
    X1=PRO(I1,2)
    X2=-24.8*X1*X1*X1*X1+31.9*X1*X1*X1-16.4*X1*X1+4.28*X1
    X3=-33.7*X2*X2*X2*X2+17.0*X2*X2*X2-4.13*X2*X2+0.327*X2+1.0
        IF(X3.LE.0.05)THEN
        X3=0.05
        END IF
    PRO(I1,2)=X2
    PRO(I1,1)=1.05*X3*PRO(I1,1)
    
    I1=7
    X1=PRO(I1,2)
    X2=-24.8*X1*X1*X1*X1+31.9*X1*X1*X1-16.4*X1*X1+4.28*X1
    X3=-33.7*X2*X2*X2*X2+17.0*X2*X2*X2-4.13*X2*X2+0.327*X2+1.0
    IF(X3.LE.0.05)THEN
    X3=0.05
    END IF
    PRO(I1,2)=X2
    PRO(I1,1)=1.05*X3*PRO(I1,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     �o�l�萔�����肷��
!     CALCULATE SPRING STIFFNESS
      DO 10 I1=1,NPHASE
!       �΂˂̐ݒ肳��Ă���ʂɑ΂��Ă̂݌v�Z
!       CALCULATION IS CONDUCTED ONLY TO SPRING SET FACE
        IF(IBANE(I1).EQ.0)THEN
	    H=HH(I1,1)+HH(I1,2)
	    M1=IELPHA(I1,1)
	    M2=IELPHA(I1,2)
!��	    
	    K1=KIND0(M1)
	    IF(K1.EQ.8)THEN
        K1=2
        END IF	    
	    K2=KIND0(M2)
	    IF(K2.EQ.8)THEN
        K2=2
        END IF
        
        IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.8.OR. &         !!8�Ԃ�1�Ԃ��ڐG����ʂ̃o�l�͗͂�`���Ȃ�(1�Ԃ����̂܂܂ő��v�H)����
         KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN
        BANE(I1,1)=1
        BANE(I1,2)=1  	    
        GOTO 10
         END IF               

            
        IF(KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.3.OR. & 
           KIND0(IELPHA(I1,1)).EQ.3.AND.KIND0(IELPHA(I1,2)).EQ.9)THEN
            E=(PRO(9,1)*HH(I1,1)+PRO(9,1)*HH(I1,2))/H;
	        P=(PRO(9,2)*HH(I1,1)+PRO(9,2)*HH(I1,2))/H;        
        
            BANE(I1,1)=E/(1.0+P)    
	        BANE(I1,2)=E*(1.0-P)/(1.0+P)/(1.0-2.0*P)        
            GOTO 10
           END IF

           
        IF(KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. & 
           KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.9)THEN
            E=(PRO(2,1)*HH(I1,1)+PRO(2,1)*HH(I1,2))/H;
	        P=(PRO(2,2)*HH(I1,1)+PRO(2,2)*HH(I1,2))/H;        
        
            BANE(I1,1)=E/(1.0+P)    
	        BANE(I1,2)=E*(1.0-P)/(1.0+P)/(1.0-2.0*P)        
            GOTO 10
           END IF
           
	    E=(PRO(K1,1)*HH(I1,1)+PRO(K2,1)*HH(I1,2))/H;
	    P=(PRO(K1,2)*HH(I1,1)+PRO(K2,2)*HH(I1,2))/H;

!       �o�l�萔�͕��ʂЂ��ݏ�Ԃ����肵�Ă���
!       PLANE STRAIN CONDITION IS ASSUMED
!           BANE(I1,1)= XY�i����f�j�o�l�萔
!           BANE(I1,2)= Z �i�����j�o�l�萔
!           BANE(I1,1)= XY  SHEAR SPRING STIFFNESS
!           BANE(I1,2)= Z   NORMAL SPRING STIFFNESS

        IF(KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.4.OR. &
           KIND0(IELPHA(I1,1)).EQ.4.AND.KIND0(IELPHA(I1,2)).EQ.9)THEN
        
            BANE(I1,1)=E/(1.0+P)
            BANE(I1,2)=0.001


            GOTO 10
        END IF

        BANE(I1,1)=E/(1.0+P)    
	    BANE(I1,2)=E*(1.0-P)/(1.0+P)/(1.0-2.0*P)

  
     
      END IF
10    CONTINUE 

      RETURN
    END SUBROUTINE SPRINGSTIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE RANDOMSPRING(MTYPE,TEN,NPHASE,IBANE,IELPHA,PRO,BANE,KIND0,XYZ)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::NPHASE,MTYPE
      INTEGER,INTENT(IN)::IBANE(:),IELPHA(:,:),KIND0(:)
      REAL(KIND=DP),INTENT(OUT):: BANE(:,:)
      REAL(KIND=DP),INTENT(IN)::XYZ(:,:),PRO(:,:)
      REAL(KIND=DP),INTENT(OUT)::TEN(:)
      
      INTEGER I1,KP,K,IS,J
      REAL(KIND=DP) M,S,XNOR,XX
                 
!     �o�l�萔�ƍő勭�x�ɂ΂����^����
!     GIVE VARIATION OF SPRING STIFFNESS AND STRENGTH
!!      OPEN(90,FILE="STIFF.TXT")

      DO I1=1,NPHASE
         !   �΂˂̐ݒ肳��Ă���ʂɑ΂��Ă̂݌v�Z
         !   CALCULATION IS CONDUCTED ONLY TO SPRING SET FACE
         IF(IBANE(I1).EQ.0)THEN
            !  �����^��-�����^���C�S��-�S�؁C�����^��-�S�؁C����-���ށC�����^��-���ށC�S��-����
            IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.6.AND.KIND0(IELPHA(I1,2)).EQ.6.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.6.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.6.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.7.AND.KIND0(IELPHA(I1,2)).EQ.7.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.7.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.7.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.6.OR.  &
                 KIND0(IELPHA(I1,2)).EQ.6.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.8.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.8.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.9.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.8.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.8.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.9.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.9.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.5.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.5.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.6.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.6.AND.KIND0(IELPHA(I1,2)).EQ.5.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.7.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.7.AND.KIND0(IELPHA(I1,2)).EQ.5)THEN
               IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.5.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.5.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN       !! MORTAL-MORTAL
                  KP=1
               ELSE IF(KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. &
                     KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.8.OR. &
                     KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.9.OR. &
                     KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.8.OR. &
                     KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. &
                     KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.9.OR. &
                     KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. &
                     KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.9.OR. &
                     KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.8)THEN  !! REBAR-REBAR
                     KP=2
               ELSE IF(KIND0(IELPHA(I1,1)).EQ.6.AND.KIND0(IELPHA(I1,2)).EQ.6.OR. &
                     KIND0(IELPHA(I1,1)).EQ.7.AND.KIND0(IELPHA(I1,2)).EQ.7)THEN  !! AGG-AGG
                     KP=6
               ELSE IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.  &
                    KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                    KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.8.OR.  &
                    KIND0(IELPHA(I1,1)).EQ.8.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                    KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.9.OR.  &
                    KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.1.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.5.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.8.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.5.AND.KIND0(IELPHA(I1,2)).EQ.9.OR.  &
                 KIND0(IELPHA(I1,1)).EQ.9.AND.KIND0(IELPHA(I1,2)).EQ.5)THEN
                    KP=10 
               ELSE 
                  KP=4
               END IF
               !      �e�����ނ̏ꍇ�͂΂������
               !      FOR ELASTIC AGGREGATE AND REBAR, NO VARIATION
               IF(MTYPE.EQ.61.OR.MTYPE.EQ.62)THEN
                  IF(KP.EQ.2)THEN
                     TEN(I1)=PRO(2,3)
                  ELSE IF(KP.EQ.6)THEN
                     TEN(I1)=PRO(6,3)
                  !!   TEN(I1)=10000000
                     GOTO 100
                  END IF
               END IF
               !     GENERATE RANDOM VALUE
               !     VARIATION IS BASED ON NORMAL DISTRIBUTION
               !     AVERAGE 0, BETWEEN 3 AND -3, 99% VALUES ARE SET IN RAND
               K=12
               M=0
               S=1.0
               XNOR=0.0
               DO J=1,K
                  XNOR=XNOR!+RAND(0)
               END DO
               XNOR=XNOR-6.0
               XX=XNOR*S+M
               !     2������͂Ɠ��l�̂΂���ɂ���
               !     SAME VARIATION AS 2D ANALYSIS
               XX=XX*(-0.2*PRO(KP,3)+1.5)
               IF(PRO(KP,3).GT.7.5)THEN
                  XX=0
               END IF
               IF(XX.LE.-2.999)THEN
                  XX=-2.999
               END IF
!         READ(IRAN1,888) XX
!         �΂��
!         VARIATION
               XX=0
               IF(KP.EQ.5)THEN
                  TEN(I1)=2*PRO(1,3)/3+PRO(KP,3)*XX/3
               ELSE
                  TEN(I1)=PRO(KP,3)+PRO(KP,3)*XX/3
               END IF
               BANE(I1,1)=(BANE(I1,1)+BANE(I1,1)*XX/3)
               BANE(I1,2)=(BANE(I1,2)+BANE(I1,2)*XX/3)
100       CONTINUE          
            ELSE
!		  TEN(I1)=10000000
            END IF
         END IF
      END DO

!!      CLOSE(90)

      RETURN
    END SUBROUTINE RANDOMSPRING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE SKY(ISKY,STDEF,STDEF1,STDEFBI,STDEF1BI,NXFIX,NYFIX,NZFIX,NSFIX,IXFIX,&
         IYFIX,IZFIX,ISFIX,NXDISP,NYDISP,NZDISP,IXDISP,IYDISP,IZDISP,IXDISPI,IYDISPI,&
         IZDISPI,XDISP,YDISP,ZDISP,NELE)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::IXFIX(:),IYFIX(:),IZFIX(:),ISFIX(:),IXDISP(:),IYDISP(:),IZDISP(:),&
           IXDISPI(:),IYDISPI(:),IZDISPI(:),NELE,NXFIX,NYFIX,NZFIX,NSFIX,NXDISP,NYDISP,NZDISP
      REAL(KIND=DP),INTENT(IN)::XDISP(:),YDISP(:),ZDISP(:)
      INTEGER,INTENT(OUT)::ISKY(:)
      REAL(KIND=DP),INTENT(OUT)::STDEF(:),STDEF1(:),STDEFBI(:),STDEF1BI(:)
      
      INTEGER I1,IA,II

!     �x�_����
!     FIXED ELEMENT CONDITION
      DO  I1=1,NXFIX
         IA=IXFIX(I1)
         ISKY((IA-1)*6+1)=1
      END DO
      DO  I1=1,NYFIX
         IA=IYFIX(I1)
         ISKY((IA-1)*6+2)=1
      END DO
      DO  I1=1,NZFIX
         IA=IZFIX(I1)
         ISKY((IA-1)*6+3)=1
      END DO
      DO  I1=1,NSFIX
         IA=ISFIX(I1)
         ISKY((IA-1)*6+4)=1
         ISKY((IA-1)*6+5)=1
         ISKY((IA-1)*6+6)=1
      END DO
      !     �����ψʏ���
      !     FORCE DISPLACEMENT ELEMENT CONDITION
      DO I1=1,NXDISP
         IA=IXDISP(I1)
         ISKY((IA-1)*6+1)=1
         STDEF((IA-1)*6+1)=XDISP(I1)
         IF(IXDISPI(I1).EQ.3.OR.IXDISPI(I1).EQ.4)THEN	  
            STDEFBI((IA-1)*6+1)=0
         ELSE
            STDEFBI((IA-1)*6+1)= XDISP(I1) !change25112015
         END IF
      END DO
      DO I1=1,NYDISP
         IA=IYDISP(I1)
         ISKY((IA-1)*6+2)=1
         STDEF((IA-1)*6+2)=YDISP(I1)
         IF(IYDISPI(I1).EQ.3.OR.IYDISPI(I1).EQ.4)THEN 
            STDEFBI((IA-1)*6+2)=0
         ELSE
            STDEFBI((IA-1)*6+2)= YDISP(I1)!change25112015
         END IF
      END DO
      DO I1=1,NZDISP
         IA=IZDISP(I1)
         ISKY((IA-1)*6+3)=1
         STDEF((IA-1)*6+3)=ZDISP(I1)
         IF(IZDISPI(I1).EQ.3.OR.IZDISPI(I1).EQ.4)THEN  
            STDEFBI((IA-1)*6+3)=0
         ELSE
            STDEFBI((IA-1)*6+3)=ZDISP(I1) !change25112015
         END IF
      END DO
      !     ISKY���v�Z�p�ɕό`
      !     MODIFICATION OF ISKY
      !     STDEF1���쐬
      !     MAKE STDEF1
      II=0
      DO I1=1,6*NELE
         IF(ISKY(I1).EQ.0)THEN
            ISKY(I1)=-II
         ELSE
            II=II+1
            ISKY(I1)=II
            STDEF1(II)=STDEF(I1)
            STDEF1BI(II)=STDEFBI(I1)
         END IF
      END DO

      RETURN
    END SUBROUTINE SKY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE KLOCATION(NELE,ITMAT,ISKY,NPHELEM,IPHELEM,IELPHA,IBANE,IK11,IK12,NK11,NK12)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::NELE,ITMAT
      INTEGER,INTENT(IN)::ISKY(:),NPHELEM(:),IPHELEM(:,:),IELPHA(:,:),IBANE(:)
      INTEGER,INTENT(OUT)::IK11(:),IK12(:),NK11(:,:),NK12(:,:)

      INTEGER::IPO1(ITMAT),IPO2(18)
      INTEGER I1,I2,I3,I4,I11,I12,IE1,IE2,IEE,IA,IA1,IA2,IB,IPH
      
      !     IK11   NUMBER OF NON-ZERO CELL ON EACH LINE OF K11 
      !     NK11   MEMORIZE COLUMN NUMBER OF NON-ZERO CELL ON EACH LINE OF K11

      DO I1=1,6*NELE
         IK11(I1)=0
         IK12(I1)=0
      END DO
      !     DIAGONAL MATRIX
      DO I1=1,NELE
         DO I2=(I1-1)*6+1,I1*6
            IF(ISKY(I2).LE.0)THEN
               I11=I2+ISKY(I2)
               DO I3=(I1-1)*6+1,I1*6
                  IF(ISKY(I3).LE.0)THEN
                     I12=I3+ISKY(I3)
                     IK11(I11)=IK11(I11)+1
                     NK11(I11,IK11(I11))=I12
                  ELSE
                     I12=ISKY(I3)
                     IK12(I11)=IK12(I11)+1
                     NK12(I11,IK12(I11))=I12
                  END IF
               END DO
            END IF
         END DO
      END DO
      !     OTHER BANDS
      DO I1=1,NELE
         DO I2=1,NPHELEM(I1)
            IPH=IPHELEM(I1,I2)  
            IF(IBANE(IPH).EQ.0)THEN
               IE1=IELPHA(IPH,1)
               IE2=IELPHA(IPH,2)
               IF(IE1.EQ.I1)THEN
                  IEE=IE2
               ELSE
                  IEE=IE1
               END IF
               DO I3=(I1-1)*6+1,I1*6
                  IF(ISKY(I3).LE.0)THEN
                     I11=I3+ISKY(I3)
                     DO I4=(IEE-1)*6+1,IEE*6
                        IF(ISKY(I4).LE.0)THEN
                           I12=I4+ISKY(I4) 
                           IK11(I11)=IK11(I11)+1
                           NK11(I11,IK11(I11))=I12
                        ELSE
                           I12=ISKY(I4)
                           IK12(I11)=IK12(I11)+1
                           NK12(I11,IK12(I11))=I12
                        END IF
                     END DO
                  END IF
               END DO
            END IF
         END DO
      END DO
      !     ARRANGE ORDER OF NK11
      IB=ABS(ISKY(6*NELE))
      DO I1=1,6*NELE-IB
         DO I2=1,ITMAT
            IPO1(I2)=0
         END DO
         IA=1
         IPO1(1)=NK11(I1,1)
         DO I2=2,IK11(I1)
            IA1=NK11(I1,I2)
            DO I3=1,IA
               IA2=IPO1(I3)
               IF(IA1.LE.IA2)THEN
                  DO I4=IA,I3,-1
                     IPO1(I4+1)=IPO1(I4)
                  END DO
                  IPO1(I3)=IA1
                  IA=IA+1
                  GOTO 110
               END IF
               IF(I3.EQ.IA)THEN
                  IPO1(I3+1)=IA1
                  IA=IA+1
                  GOTO 110
               END IF
            END DO
110      END DO
         !       COMPLETE THE ARRANGEMENT UNTIL I1 LINE, RE-RECORD
         DO I2=1,IK11(I1)
            NK11(I1,I2)=IPO1(I2)
         END DO
      END DO
      !     ARRANGE THE ORDER OF NK12
      DO I1=1,6*NELE-IB
         DO I2=1,18
            IPO2(I2)=0
         END DO
         IA=1
         IPO2(1)=NK12(I1,1)
         DO 170 I2=2,IK12(I1)
            IA1=NK12(I1,I2)
            DO I3=1,IA
               IA2=IPO2(I3)
               IF(IA1.LE.IA2)THEN
                  DO I4=IA,I3,-1
                     IPO2(I4+1)=IPO2(I4)
                  END DO
                  IPO2(I3)=IA1
                  IA=IA+1
                  GOTO 170
               END IF
               IF(I3.EQ.IA)THEN
                  IPO2(I3+1)=IA1
                  IA=IA+1
                  GOTO 170
               END IF
            END DO
170      END DO
         !       COMPLETE THE ARRANGEMENT UNTIL I1 LINE, RE-RECORD
         DO I2=1,IK12(I1)
            NK12(I1,I2)=IPO2(I2)
         END DO
      END DO

      RETURN
    END SUBROUTINE KLOCATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MAKEMATRIX(NPHASE,NELE,COD,AREA,BANE,IBANE,NNOPHA,INOPHA,&
         IELPHA,ISKY,NK11,XK11,NK12,XK12,IK11,IK12,XYZ,HH,CPHASE,ITMAT)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::NPHASE,NELE,ITMAT
      INTEGER,INTENT(IN)::IBANE(:),NNOPHA(:),INOPHA(:,:),IELPHA(:,:),ISKY(:)
      REAL(KIND=DP),INTENT(IN)::COD(:,:),AREA(:),XYZ(:,:),BANE(:,:),CPHASE(:,:),HH(:,:)
      REAL(KIND=DP),INTENT(OUT)::XK11(:,:),XK12(:,:)
      INTEGER,INTENT(IN)::NK11(:,:),NK12(:,:),IK11(:),IK12(:)
     
      INTEGER I1,I2,IS
      REAL(KIND=DP)::EMAT(12,12)
      INTEGER,ALLOCATABLE::ICK1(:,:),ICK2(:,:)

      ALLOCATE(ICK1(6*NELE,ITMAT))
      ALLOCATE(ICK2(6*NELE,18))
      XK11 =0.D0
      ICK1 =0
      XK12 =0.D0
      ICK2 =0
      ! MAKE GLOBAL MATRIX
      DO IS=1,NPHASE
         IF(IBANE(IS).EQ.0)THEN
            ! MAKE LOCAL MATRIX
            CALL ELEMATRIX(COD,AREA,BANE,NNOPHA,INOPHA,IELPHA,EMAT,XYZ,HH,&
                 CPHASE,IS)
            !       INSTALL LOCAL MATRIX TO GLOBAL MATRIX
            CALL MKTMAT(NELE,IELPHA,ISKY,NK11,XK11,NK12,XK12,&
                 IK11,IK12,EMAT,ITMAT,IS,ICK1,ICK2)
         END IF
      END DO

!     DELETE DYNAMIC ARRAY
      DEALLOCATE (ICK1)
      DEALLOCATE (ICK2)

      RETURN
    END SUBROUTINE MAKEMATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE ELEMATRIX(COD,AREA,BANE,NNOPHA,INOPHA,IELPHA,EMAT,XYZ,HH,CPHASE,IS)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::IS
      INTEGER,INTENT(IN)::NNOPHA(:),INOPHA(:,:),IELPHA(:,:)
      REAL(KIND=DP),INTENT(IN)::COD(:,:),AREA(:),CPHASE(:,:),BANE(:,:),XYZ(:,:),HH(:,:)
      REAL(KIND=DP),INTENT(OUT)::EMAT(:,:)
      
      REAL(KIND=DP)::BMAT1(3,12)

!     CALCULATION OF B MATRIX
      CALL BMATRIX0(COD,NNOPHA,INOPHA,IELPHA,CPHASE,XYZ,BMAT1,IS)
!     MAKE LOCAL MATRIX
      CALL MKELEMAT(AREA,HH,BANE,EMAT,BMAT1,IS)

      RETURN
    END SUBROUTINE ELEMATRIX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE BMATRIX0(COD,NNOPHA,INOPHA,IELPHA,CPHASE,XYZ,BMAT1,IS)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::IS
      INTEGER,INTENT(IN)::NNOPHA(:),INOPHA(:,:),IELPHA(:,:)
      REAL(KIND=DP),INTENT(IN)::COD(:,:),CPHASE(:,:),XYZ(:,:)
      REAL(KIND=DP),INTENT(OUT)::BMAT1(:,:)
      
      INTEGER IA,IB1,IB2,IFLAG,I1
      REAL(KIND=DP) XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,X1,Y1,Z1,X2,Y2,Z2,A,B,C,D,& 
           HH2,XG1,YG1,ZG1,AUX,AUY,AUZ,AVX,AVY,AVZ,ACX,ACY,ACZ,AWX,AWY,AWZ,&
           XLU,XLV,XLW,EUX,EUY,EUZ,EVX,EVY,EVZ,EWX,EWY,EWZ,XG2,YG2,ZG2,FLAG,XX,YY,ZZ
      !     CALCULATION OF B MATRIX
      !     CALCULATE NORMAL VECTOR USING 
      !        THE TRIANGLE CONSISTS OF LINE OF FACE AND CENTER OF GRAVITY OF FACE
      IA=NNOPHA(IS)
      XA=CPHASE(IS,1)
      YA=CPHASE(IS,2)
      ZA=CPHASE(IS,3)
      DO I1=1,IA
         IF(I1.NE.IA)THEN
            IB1=INOPHA(IS,I1)
            XB=COD(IB1,1)
            YB=COD(IB1,2)
            ZB=COD(IB1,3)
            IB2=INOPHA(IS,I1+1)
            XC=COD(IB2,1)
            YC=COD(IB2,2)
            ZC=COD(IB2,3)
         ELSE
            IB1=INOPHA(IS,I1)
            XB=COD(IB1,1)
            YB=COD(IB1,2)
            ZB=COD(IB1,3)
            IB2=INOPHA(IS,1)
            XC=COD(IB2,1)
            YC=COD(IB2,2)
            ZC=COD(IB2,3)
         END IF
!        �x�N�g��
         X1=-(XA-XB)
         Y1=-(YA-YB)
         Z1=-(ZA-ZB)
         X2=-(XA-XC)
         Y2=-(YA-YC)
         Z2=-(ZA-ZC)
!        NORMAL VECTOR
         A=Y1*Z2-Z1*Y2
         B=Z1*X2-X1*Z2
         C=X1*Y2-Y1*X2
!        IF HH2 IS NOT ZERO, IT BECOMES TWO VECTORS WITH AREA
!        USE THE COORDINATE, MAKE B MATRIX
         HH2=A*A+B*B+C*C
         IF(HH2.GT.0)THEN
            GOTO 300
         END IF
      END DO

300   CONTINUE
!     FLAG FOR DECISION OF DIRECTION OF W VECTOR
!     SET THE DIRECTION TO MAKE ELEMENT 1 BE UNDER THE W VECTOR DIRECTION 

!     AX+BY+CZ+D=0�Ƃ����
!     IF AX+BY+CZ+D=0
      D=-A*XA-B*YA-C*ZA
!     CENTER OF GRAVITY OF ELEMENT 1
      IA=IELPHA(IS,1)
      XG1=XYZ(IA,1)
      YG1=XYZ(IA,2)
      ZG1=XYZ(IA,3)
      FLAG=A*XG1+B*YG1+C*ZG1+D
!     IF FLAG<0, W DOES NOT CHANGE
!     IF FLAG>0, W*-1
      IF(FLAG.LT.0)THEN
         IFLAG=0
      ELSE 
         IFLAG=-1
      END IF
!     CALCULATION OF R MATRIX
      AUX=XB-XA
      AUY=YB-YA
      AUZ=ZB-ZA
      ACX=XC-XA
      ACY=YC-YA
      ACZ=ZC-ZA
      AWX=AUY*ACZ-AUZ*ACY
      AWY=AUZ*ACX-AUX*ACZ
      AWZ=AUX*ACY-AUY*ACX
!     CANGHE THE DIRECITON OF W DEPENDS ON IFLAG
      IF(IFLAG.EQ.-1)THEN
         AWX=-AWX
         AWY=-AWY
         AWZ=-AWZ
      END IF
      AVX=AWY*AUZ-AWZ*AUY
      AVY=AWZ*AUX-AWX*AUZ
      AVZ=AWX*AUY-AWY*AUX
      XLU=SQRT(AUX**2+AUY**2+AUZ**2)
      XLV=SQRT(AVX**2+AVY**2+AVZ**2)
      XLW=SQRT(AWX**2+AWY**2+AWZ**2)
      EUX=AUX/XLU
      EUY=AUY/XLU
      EUZ=AUZ/XLU
      EVX=AVX/XLV
      EVY=AVY/XLV
      EVZ=AVZ/XLV
      EWX=AWX/XLW
      EWY=AWY/XLW
      EWZ=AWZ/XLW
!     CALCULATION OF B MATRIX
!     CENTER OF COORDINATE OF FACE
      XX=CPHASE(IS,1)
      YY=CPHASE(IS,2)
      ZZ=CPHASE(IS,3)
!     CENTER OF COORDINATE OF ELEMENT 1
      IA=IELPHA(IS,1)
      XG1=XYZ(IA,1)
      YG1=XYZ(IA,2)
      ZG1=XYZ(IA,3)
!     CENTER OF COORDINATE OF ELEMENT 2
      IA=IELPHA(IS,2)
      XG2=XYZ(IA,1)
      YG2=XYZ(IA,2)
      ZG2=XYZ(IA,3)
!     CALCULATION OF B MATRIX
      BMAT1(1,1)=-EUX
      BMAT1(2,1)=-EVX
      BMAT1(3,1)=-EWX
      BMAT1(1,2)=-EUY
      BMAT1(2,2)=-EVY
      BMAT1(3,2)=-EWY
      BMAT1(1,3)=-EUZ
      BMAT1(2,3)=-EVZ
      BMAT1(3,3)=-EWZ
      BMAT1(1,4)=-(EUY*(ZZ-ZG1)-EUZ*(YY-YG1))
      BMAT1(2,4)=-(EVY*(ZZ-ZG1)-EVZ*(YY-YG1))
      BMAT1(3,4)=-(EWY*(ZZ-ZG1)-EWZ*(YY-YG1))
      BMAT1(1,5)=-(EUZ*(XX-XG1)-EUX*(ZZ-ZG1))
      BMAT1(2,5)=-(EVZ*(XX-XG1)-EVX*(ZZ-ZG1))
      BMAT1(3,5)=-(EWZ*(XX-XG1)-EWX*(ZZ-ZG1))
      BMAT1(1,6)=-(EUX*(YY-YG1)-EUY*(XX-XG1))
      BMAT1(2,6)=-(EVX*(YY-YG1)-EVY*(XX-XG1))
      BMAT1(3,6)=-(EWX*(YY-YG1)-EWY*(XX-XG1))      

      BMAT1(1,7)=EUX
      BMAT1(2,7)=EVX
      BMAT1(3,7)=EWX
      BMAT1(1,8)=EUY
      BMAT1(2,8)=EVY
      BMAT1(3,8)=EWY
      BMAT1(1,9)=EUZ
      BMAT1(2,9)=EVZ
      BMAT1(3,9)=EWZ
      BMAT1(1,10)=(EUY*(ZZ-ZG2)-EUZ*(YY-YG2))
      BMAT1(2,10)=(EVY*(ZZ-ZG2)-EVZ*(YY-YG2))
      BMAT1(3,10)=(EWY*(ZZ-ZG2)-EWZ*(YY-YG2))
      BMAT1(1,11)=(EUZ*(XX-XG2)-EUX*(ZZ-ZG2))
      BMAT1(2,11)=(EVZ*(XX-XG2)-EVX*(ZZ-ZG2))
      BMAT1(3,11)=(EWZ*(XX-XG2)-EWX*(ZZ-ZG2))
      BMAT1(1,12)=(EUX*(YY-YG2)-EUY*(XX-XG2))
      BMAT1(2,12)=(EVX*(YY-YG2)-EVY*(XX-XG2))
      BMAT1(3,12)=(EWX*(YY-YG2)-EWY*(XX-XG2))
      
      RETURN
    END SUBROUTINE BMATRIX0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MKELEMAT(AREA,HH,BANE,EMAT,BMAT1,IS)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::IS
      REAL(KIND=DP),INTENT(IN)::AREA(:),HH(:,:),BANE(:,:),BMAT1(:,:)
      REAL(KIND=DP),INTENT(OUT)::EMAT(:,:)
    
      REAL(KIND=DP)::BMAT2(3,12)
      REAL(KIND=DP) H,DD,ARE,XKX,XKY,XKZ
      INTEGER I1,I2,I3
      
!     K_ELEM=[B]t[DA][B]
      ARE=AREA(IS)
      H=HH(IS,1)+HH(IS,2)
      XKX=BANE(IS,1)*ARE/H
      XKY=BANE(IS,1)*ARE/H
      XKZ=BANE(IS,2)*ARE/H
      DO I1=1,12
         BMAT2(1,I1)=BMAT1(1,I1)*XKX
         BMAT2(2,I1)=BMAT1(2,I1)*XKY
         BMAT2(3,I1)=BMAT1(3,I1)*XKZ
      END DO
!     MATRIX SIZE IS 12 X 12
      DO I1=1,12
         DO I2=1,12
            DD=0.0
            DO I3=1,3
               DD=DD+BMAT2(I3,I1)*BMAT1(I3,I2)
            END DO
            EMAT(I1,I2)=DD
         END DO
      END DO
      
      RETURN
    END SUBROUTINE MKELEMAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE MKTMAT(NELE,IELPHA,ISKY,NK11,XK11,NK12,XK12,IK11,IK12,EMAT,ITMAT,IS,ICK1,ICK2)
      
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::ITMAT,IS,NELE
      INTEGER,INTENT(IN)::IELPHA(:,:),ISKY(:),NK11(:,:),NK12(:,:),IK11(:),IK12(:)
      REAL(KIND=DP),INTENT(IN)::EMAT(:,:)
      REAL(KIND=DP),INTENT(OUT)::XK11(:,:),XK12(:,:)
      INTEGER,INTENT(OUT)::ICK1(:,:),ICK2(:,:)
      INTEGER::IE(2)
      REAL(KIND=DP)::EEMAT(6,6)
      
      INTEGER I1,I2,I3,I4,I5,IA1,EX,EY,I11,I12

!     INSTALL LOCAL MATRIX TO XK11 AND XK12
!     ACCORD THE LOCATION TO NK11 AND NK12
!     IS IS FACE NUMBER
      IE(1)=IELPHA(IS,1)
      IE(2)=IELPHA(IS,2)

      DO I1=1,2
         DO I2=1,2
!        MODIFIED STRAGE WAY
            IF(I1.EQ.1.AND.I2.EQ.1)THEN
               DO I3=1,6
                  DO I4=1,I3
                     EEMAT(I3,I4)=EMAT(I3,I4)
                  END DO
               END DO
               DO I3=1,5
                  DO I4=I3+1,6
                     EEMAT(I3,I4)=EEMAT(I4,I3)
                  END DO
               END DO
            ELSE IF(I1.EQ.1.AND.I2.EQ.2)THEN
               DO I3=1,6
                  DO I4=1,6
                     EEMAT(I4,I3)=EMAT(6+I3,I4)
                  END DO
               END DO
            ELSE IF(I1.EQ.2.AND.I2.EQ.1)THEN
               DO I3=1,6
                  DO I4=1,6
                     EEMAT(I3,I4)=EMAT(6+I3,I4)
                  END DO
               END DO
            ELSE
               DO I3=1,6
                  DO I4=1,I3
                     EEMAT(I3,I4)=EMAT(6+I3,6+I4)
                  END DO
               END DO
               DO I3=1,5
                  DO I4=I3+1,6
                     EEMAT(I3,I4)=EEMAT(I4,I3)
                  END DO
               END DO
            END IF
            EX=(IE(I1)-1)*6
            EY=(IE(I2)-1)*6
            DO I3=1,6
               IF(ISKY(EX+I3).LE.0)THEN
                  I11=EX+I3+ISKY(EX+I3)
                  DO I4=1,6
                     IF(ISKY(EY+I4).LE.0)THEN
                        I12=EY+I4+ISKY(EY+I4)
                        !   SEARCH (III,?)=I12 OF NK11
                        DO I5=1,IK11(I11)
                           IA1=NK11(I11,I5)
                           IF(IA1.EQ.I12)THEN
                              XK11(I11,I5)=XK11(I11,I5)+EEMAT(I3,I4)
                              ICK1(I11,I5)=ICK1(I11,I5)+1
                              GOTO 80
                           END IF
                        END DO
                     ELSE
                        I12=ISKY(EY+I4)
                        !   SEARCH (III,?)=I12 OF NK12
                        DO I5=1,IK12(I11)
                           IA1=NK12(I11,I5)
                           IF(IA1.EQ.I12)THEN
                              XK12(I11,I5)=EEMAT(I3,I4)
                              ICK2(I11,I5)=ICK2(I11,I5)+1
                              GOTO 80
                           END IF
                        END DO
                     END IF
80                END DO
               END IF
            END DO
         END DO
      END DO

      RETURN
    END SUBROUTINE MKTMAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE CALNZND(NNZ,NELE,NK11,XK11,IK11,ITMAT,NNN,ISKY)
      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::NELE,ITMAT
      INTEGER,INTENT(IN)::ISKY(:),NK11(:,:),IK11(:)
      REAL(KIND=DP),INTENT(IN)::XK11(:,:)
      INTEGER,INTENT(OUT)::NNZ,NNN
      
      INTEGER I1,I2,II

!     CALCULATE NNZ AND NNN FOR SOLVER
!     CALCULATE NNN
      NNN=6*NELE-ABS(ISKY(6*NELE))
!     CALCULATE NNZ
      NNZ=0
      DO I1=1,6*NELE
         IF(IK11(I1)/=0)THEN
            II=I1+ISKY(I1)
            DO I2=1,IK11(I1)
               IF(NK11(I1,I2) >= II)THEN
                  NNZ=NNZ+1
               END IF
            END DO
         ENDIF
      END DO

      RETURN
    END SUBROUTINE CALNZND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
SUBROUTINE SIMULATION(NNODE,NELE,NPHASE,COD,NNOPHA,INOPHA,NPHELEM,IPHELEM,IELPHA,KIND0,&
     NXDISP,NYDISP,NZDISP,IXDISP,IYDISP,IZDISP,IXDISPI,IYDISPI,IZDISPI,XDISP,YDISP,ZDISP,&
     IFSTEP,XYZ,AREA,HH,IBANE,CPHASE,MTYPE,IBTYPE,PRO,ITMAT,BANE,ISKY,STDEF,STDEF1,STDEFBI,&
     STDEF1BI,IK11,IK12,NNN,TEN,NNZ,ISTEPBI,NK11,NK12,XK11,XK12)
  IMPLICIT NONE
      
  INTEGER,INTENT(IN)::NNODE,NELE,NPHASE,IFSTEP,MTYPE,NXDISP,NYDISP,NZDISP,& 
       IBTYPE,ITMAT,ISTEPBI,NNN,NNZ
  INTEGER,INTENT(IN)::NNOPHA(:),INOPHA(:,:),NPHELEM(:),IPHELEM(:,:),IELPHA(:,:),KIND0(:),&
       IXDISP(:),IYDISP(:),IZDISP(:),IXDISPI(:),IYDISPI(:),IZDISPI(:),IBANE(:),ISKY(:)
  INTEGER,INTENT(IN)::NK11(:,:),NK12(:,:),IK11(:),IK12(:)
  REAL(KIND=DP),INTENT(IN)::COD(:,:),AREA(:),HH(:,:),XYZ(:,:),XDISP(:),YDISP(:),ZDISP(:),&
       CPHASE(:,:),PRO(:,:),TEN(:),BANE(:,:),STDEF(:),STDEFBI(:)
  REAL(KIND=DP),INTENT(IN)::STDEF1(:),STDEF1BI(:),XK11(:,:),XK12(:,:)
  REAL(KIND=DP),ALLOCATABLE::FA(:),U(:),DU1(:),DU(:),DDF(:),STRA(:,:),STRES(:,:), &
       HIS(:,:),HIST(:,:),SHIS(:),SHIST(:),HISP(:,:),HISN(:,:),HISTP(:,:),HISTN(:,:)
  INTEGER,ALLOCATABLE::JYD(:,:),JCRA(:,:)   
  INTEGER ISTEP,FLAG,ITER,I,AA,I1
  INTEGER IOUT1,IOUT3,IOUT4,IOUT5,IOUT6,IOUT7,IOUT11,IOUT12,IOUT13,IOUT14,IOUT15,IOUT16,&
       IOUT17,IOUT18,IOUT19,IOUT20,IOUT21,IOUT22,IOUT23,IOUT24,IOUT25,IOUT26,IOUT27,IREC,IREC1
  REAL(KIND=DP) RN,SN,RATE,JV,STEPFLAG
  REAL(KIND=DP),ALLOCATABLE::a(:),BMAT(:,:,:),CWORK(:,:,:),B1(:),b(:),x(:)
  INTEGER,ALLOCATABLE::ia(:),ja(:)
  INTEGER::phase,maxfct,error,mnum,motype,nrhs,msglvl,solver,error1
  TYPE(MKL_PARDISO_HANDLE),ALLOCATABLE::pt(:)
  INTEGER,ALLOCATABLE::IPARM(:)
  INTEGER::idum(1)
  REAL(KIND=DP)::ddum(1)

  ALLOCATE(IPARM(64),pt(64))
  IPARM = 0
  DO I1=1,64
     pt(I1)%DUMMY=0
  END DO

  IPARM(1) = 1 ! no solver default
  IPARM(2) = 2 ! fill-in reordering from METIS
  IPARM(3) = 6 ! number of processors
  IPARM(4) = 0 ! no iterative-direct algorithm
  IPARM(5) = 0 ! no user fill-in reducing permutation
  IPARM(6) = 0 ! =0 solution on the first n compoments of x
  IPARM(8) = 8 ! numbers of iterative refinement steps
  IPARM(10) = 13 ! perturbe the pivot elements with 1E-13
  IPARM(11) = 1 ! use nonsymmetric permutation and scaling MPS
  IPARM(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric).
  ! Try iparm(13) = 1 in case of inappropriate accuracy
  IPARM(14) = 0 ! Output: number of perturbed pivots
  IPARM(18) = -1 ! Output: number of nonzeros in the factor LU
  IPARM(19) = -1 ! Output: Mflops for LU factorization
  IPARM(20) = 0 ! Output: Numbers of CG Iterations
  IPARM(60) = 0 ! in core memory 
  nrhs   = 1
  maxfct = 1
  mnum   = 1 
  error  = 0 ! initialize error flag
  msglvl = 0 ! print statistical information
  motype  = -2 ! symmetric, indefinite
  solver = 0 ! direct solver

!     ��͉ߒ��ۑ��@�t�@�C���쐬
!     OUTPUT FILE
  IOUT1=810
  IOUT3=910	  
  IOUT4=890
  IOUT5=870
  IOUT6=880
  IOUT11=900
  IOUT7=920
  IOUT12=930
  IOUT13=940
  IOUT14=950
  IOUT15=960
  IOUT16=970
  IOUT17=980
  IOUT18=981
  IOUT19=982
  IOUT20=983
  IOUT21=984
  IOUT22=901
  IOUT23=902
  IOUT24=903
  IOUT25=904
  IOUT26=905
  IOUT27=906

!     OPEN(800,FILE='CONVERGENCE.TXT')
!     OPEN(IOUT1,FILE='DISP-LOAD(CSV).CSV')
  OPEN(IOUT11,FILE='DISP-LOAD(TEXT).TXT')
  OPEN(IOUT22,FILE='DISP-LOAD(CSV)-RIGHT.CSV')
  OPEN(IOUT23,FILE='DISP-LOAD(TEXT)-RIGHT.TXT')
  OPEN(IOUT24,FILE='DISP-LOAD(CSV)-LEFT.CSV')
  OPEN(IOUT25,FILE='DISP-LOAD(TEXT)-LEFT.TXT')
!     OPEN(IOUT3,FILE='SOLVE2-INFO.TXT')
!     OPEN(IOUT4,FILE='UNBALANCE.TXT')
!     OPEN(IOUT5,FILE='STRA-STRES.TXT')
!     OPEN(IOUT7,FILE='TIME.TXT')
!     OPEN(IOUT6,FILE='OUTPUT.INFO',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)
  OPEN(IOUT6,FILE='OUTPUT.INFO',STATUS='UNKNOWN',FORM='BINARY')
  OPEN(IOUT12,FILE='N-BREAKSP.CSV')
  OPEN(IOUT18,FILE='CRACK0000mm.TXT')
  OPEN(IOUT13,FILE='CRACK0002mm.TXT')
  OPEN(IOUT14,FILE='CRACK0005mm.TXT')
  OPEN(IOUT15,FILE='CRACK0010mm.TXT')
  OPEN(IOUT16,FILE='CRACK0030mm.TXT')
  OPEN(IOUT19,FILE='CRACK0100mm.TXT')
  OPEN(IOUT20,FILE='CRACK0300mm.TXT')
  OPEN(IOUT21,FILE='CRACK1000mm.TXT')
  OPEN(IOUT26,FILE='CRACK00005mm.TXT')
  OPEN(IOUT27,FILE='CRACK0001mm.TXT')
  OPEN(IOUT17,FILE='AGGREGATE-SURFACE.INFO',STATUS='UNKNOWN',FORM='BINARY')
!     OPEN(IOUT17,FILE='AGGREGATE-SURFACE.INFO',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)
!     SIMULATION START TIME
  print *, "     SUBROUTINE SIMULATION"
  ALLOCATE(BMAT(NPHASE,3,12),CWORK(NPHASE,5,3))
!     MAKE B MATRIX 
  CALL BMATRIX(NPHASE,COD,NNOPHA,INOPHA,IELPHA,CPHASE,XYZ,BMAT,IBANE,CWORK)
  ALLOCATE(a(NNZ),ja(NNZ),ia(NNN+1),b(NNN),x(NNN))
!     MAKE ARRAY FOR SOLVER 
  CALL MKARRAY(NNN,IK11,ia,ja,a,XK11,NK11)
  ISTEP = 0
  FLAG = 0
  IREC = 0
  IREC1= 0
  STEPFLAG = 0
  ALLOCATE(B1(NNN))
  ALLOCATE(U(6*NELE),DU1(6*NELE),HIS(NPHASE,2),HIST(NPHASE,2),SHIS(NPHASE),SHIST(NPHASE),HISP(NPHASE,2),HISN(NPHASE,2),HISTP(NPHASE,2),HISTN(NPHASE,2))
  ALLOCATE(FA(6*NELE),DU(6*NELE),DDF(6*NELE),STRA(NPHASE,3),STRES(NPHASE,3))
  ALLOCATE(JYD(NPHASE,2),JCRA(NPHASE,10))

  DO WHILE(STEPFLAG==0)      
!!!!!!START STEP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ISTEP = ISTEP+1
     ITER = 1
     !   IN CASE STEP 1 OT THE STEP IN WHICH LOADING CONDITION IS CHANGED
     IF(ISTEP.EQ.1.OR.ISTEP.EQ.ISTEPBI+1)THEN
        !       MAKE F IN F=KU FOR ITERATION=1
        CALL FORCE1(NELE,IK12,NK12,XK12,STDEF1,STDEF1BI,B1,NNN,U,ISTEPBI,ISTEP)
        b = 0.d0
        x = 0.d0
        b(1:NNN)=B1(1:NNN)
        write(*,*)"    START SOLVE1"
        phase=11
        CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,a,ia,ja,&
                               idum,nrhs,IPARM,msglvl,ddum,ddum,error)
        phase=22
        CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,a,ia,ja,&
                               idum,nrhs,IPARM,msglvl,ddum,ddum,error)
        phase=33
        iparm(8)=0
        CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,a,ia,ja,&
                               idum,nrhs,IPARM,msglvl,b,x,error)
        phase=-1
        CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,a,ia,ja,&
                               idum,nrhs,IPARM,msglvl,ddum,ddum,error1)
        write(*,*)"    FINISHI SOLVE1"
        AA = 0
        DO I1=1,6*NELE
           IF(ISKY(I1).LE.0)THEN
              AA = AA+1
              DU1(I1) = x(AA)
           ELSE
              IF(ISTEP==1)THEN
                 DU1(I1)=STDEF(I1)
              ELSE
                 DU1(I1)=STDEFBI(I1)
              END IF
           END IF
        END DO
        
!place new du here !change22052016
!OPEN(600,FILE='DUTranAF_36.TXT',status='old',ACCESS='sequential',FORM='formatted',ACTION='read')	     
!RDU=0.D0
!DO I1=1,16146*6
!READ(600,110) RDU(I1)
!!U(I1)=RDU(I1)
!END DO
!110 format(F15.10)
!close(600)
!print *, "RDU(1)=",RDU(1)
!print *, "RDU(2)=",RDU(2)
!print *, "RDU(3)=",RDU(3)
!print *, "RDU(4)=",RDU(4)
!    
!pause



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
     !  print *, "ISTEP=1"
        !change25052016
     !   OPEN(10014,FILE='CheckDUTran23.TXT')
     !   WRITE (10014,'(F15.10)')RDU
     !   CLOSE(10014)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        
        IF(ISTEP==1)THEN
           DO I1 = 1,NPHASE
              DO I = 1,2
                 JYD(I1,I)=0
              END DO
           END DO
        END IF
    ENDIF
    !   RENEW STRESS HISTORY
    DO I = 1,NPHASE
       HIS(I,1) = HIST(I,1)
       HIS(I,2) = HIST(I,2)
       SHIS(I) = SHIST(I)
       HISN(I,1) = HISTN(I,1)
       HISN(I,2) = HISTN(I,2)
       HISP(I,1) = HISTP(I,1)       
       HISP(I,2) = HISTP(I,2)
       
    END DO
    !   CALCULATION OF STRAIN, STRESS AND UNBALANCED FORCE
    CALL ACTUALFORCE1(NELE,NPHASE,KIND0,IELPHA,FA,DDF,U,DU1,IBANE,ISKY,AREA,HH,BANE,STRES,STRA,&
         PRO,BMAT,CWORK,RN,SN,IOUT4,IOUT5,TEN,MTYPE,JYD,ISTEP,JCRA,HIS,HIST,SHIS,SHIST,XYZ,HISP,HISN,HISTP,HISTN)
    !   JUDGEMENT (JV=JUDGEMENT VALUE)
    JV=0.00001
    RATE=SN/RN
    WRITE(*,'(A,I3,A,I3)')'STEP=',ISTEP,'   ITER=',ITER
    WRITE(*,'(A,D15.5,A,D15.5,A,D15.5)')'SN=',SN,'  RN=',RN,'   RATE=',RATE
   IF(ITER.GE.1)THEN
       FLAG = 1
    ELSE
       FLAG = 0
       phase=11
      IPARM(8) = 8
      CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,a,ia,ja,&
            idum,nrhs,IPARM,msglvl,ddum,ddum,error)
       phase=22
       CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,a,ia,ja,&
            idum,nrhs,IPARM,msglvl,ddum,ddum,error)
    ENDIF
!!!!!!START ITERATION CALCULATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO WHILE(FLAG == 0)
       ITER = ITER+1
       AA = 0
       DO I1 = 1,6*NELE
          IF(ISKY(I1).LE.0)THEN
             AA = AA+1
             b(AA) = DDF(I1)
          END IF
       END DO
       x = 0.d0
       phase = 33
       iparm(8)=0
       CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,a,ia,ja,&
            idum,nrhs,IPARM,msglvl,b,x,error)
       AA = 0
       DO I1=1,6*NELE
          IF(ISKY(I1).LE.0)THEN
             AA=AA+1
             DU(I1)=x(AA)
          ELSE
             DU(I1)=0.d0
          END IF
       END DO
       CALL ACTUALFORCE2(NELE,NPHASE,KIND0,IELPHA,FA,DDF,U,DU,IBANE,&
            ISKY,AREA,HH,BANE,STRES,STRA,PRO,BMAT,CWORK,RN,SN,&
            TEN,MTYPE,JYD,ISTEP,JCRA,HIS,HIST,SHIS,SHIST,XYZ,HISP,HISN,HISTP,HISTN)
       ! CONVERGENCE JUDGEMENT
       RATE = SN/RN
       WRITE(*,'(A,I3,A,I3)')'STEP=',ISTEP,'   ITER=',ITER
       WRITE(*,'(A,D15.5,A,D15.5,A,D15.5)')'SN=',SN,'  RN=',RN,'   RATE=',RATE
       
       IF(ITER.GE.1)THEN
          FLAG = 1
          phase = -1
          CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,ddum,idum,idum,&
               idum,nrhs,IPARM,msglvl,ddum,ddum,error1)
       ELSE
       !!!!!!!!!GIVE UP CONVERGENCE!!!!!!!!!!!!!!!!
          IF(ITER.GE.1)THEN
             FLAG = 1
             phase = -1
             CALL pardiso(pt,maxfct,mnum,motype,phase,NNN,ddum,idum,idum,&
                  idum,nrhs,IPARM,msglvl,ddum,ddum,error1)
          ENDIF
       ENDIF
    END DO
!!!!!!END ITERATION CALCULATION!!!!!!!!!!!!!!!!!!!!!!!

   
!    CALL OUTPUT_LOAD()
    CALL OUTPUT(NELE,ISTEP,NXDISP,NYDISP,NZDISP,IXDISP,IYDISP,IZDISP,IXDISPI,IYDISPI,IZDISPI,&
         KIND0,IBANE,NNOPHA,INOPHA,NPHELEM,IPHELEM,IELPHA,JCRA,JYD,FA,U,XDISP,YDISP,ZDISP,COD,&
         XYZ,STRA,STRES,IOUT6,IOUT11,IOUT12,IOUT13,IOUT14,IOUT15,IOUT16,IOUT17,IOUT18,IOUT19,IOUT20,&
         IOUT21,IOUT22,IOUT23,IOUT24,IOUT25,IOUT26,IOUT27,IREC,IREC1)

    IF(ISTEP.GE.4000)THEN        ! IF(ISTEP.GE.IFSTEP)!  for step 56 crack width 40um IF(ISTEP.GE.IFSTEP)THEN 
       STEPFLAG = 1
    ELSE
       STEPFLAG = 0
    ENDIF
 END DO
   
!	  CLOSE(800)
!	  CLOSE(IOUT1)
!	  CLOSE(IOUT3)
!	  CLOSE(IOUT4)
!	  CLOSE(IOUT5)
 CLOSE(IOUT6)
!	  CLOSE(IOUT7)
 CLOSE(IOUT11)
 CLOSE(IOUT12)
 CLOSE(IOUT13)
 CLOSE(IOUT14)
 CLOSE(IOUT15)
 CLOSE(IOUT16)
 CLOSE(IOUT17)
 CLOSE(IOUT18)
 CLOSE(IOUT19)
 CLOSE(IOUT20)
 CLOSE(IOUT21)
 CLOSE(IOUT22)
 CLOSE(IOUT23)
 CLOSE(IOUT24)
 CLOSE(IOUT25)
 CLOSE(IOUT26)
 CLOSE(IOUT27)

 RETURN
END SUBROUTINE SIMULATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BMATRIX(NPHASE,COD,NNOPHA,INOPHA,IELPHA,CPHASE,XYZ,&
     BMAT,IBANE,CWORK)
  IMPLICIT NONE
      
  INTEGER,INTENT(IN)::NPHASE
  INTEGER,INTENT(IN)::NNOPHA(:),INOPHA(:,:),IELPHA(:,:),IBANE(:)
  REAL(KIND=DP),INTENT(IN)::COD(:,:),CPHASE(:,:),XYZ(:,:)
  REAL(KIND=DP),INTENT(OUT)::CWORK(:,:,:),BMAT(:,:,:)
    
  INTEGER IS,I1,IA,IB1,IB2,IFLAG
  REAL(KIND=DP) XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,X1,X2,Y1,Y2,Z1,Z2,A,B,C,D,HH2,XG1,YG1,ZG1,&
       XG2,YG2,ZG2,XX,YY,ZZ,AUX,AUY,AUZ,ACX,ACY,ACZ,AWX,AWY,AWZ,AVX,AVY,AVZ,&
       FLAG,XLU,XLW,XLV,EUX,EUY,EUZ,EVX,EVY,EVZ,EWX,EWY,EWZ
  
  !     CALCULATION OF B MATRIX
  DO IS=1,NPHASE
     IF(IBANE(IS).EQ.0)THEN     
        !    CALCULATE NORMAL VECTOR USING 
        !    THE TRIANGLE CONSISTS OF LINE OF FACE AND CENTER OF GRAVITY OF FACE
        IA=NNOPHA(IS)
        XA=CPHASE(IS,1)
        YA=CPHASE(IS,2)
        ZA=CPHASE(IS,3)
        DO I1=1,IA
           IF(I1.NE.IA)THEN
              IB1=INOPHA(IS,I1)
              XB=COD(IB1,1)
              YB=COD(IB1,2)
              ZB=COD(IB1,3)
              IB2=INOPHA(IS,I1+1)
              XC=COD(IB2,1)
              YC=COD(IB2,2)
              ZC=COD(IB2,3)
           ELSE
              IB1=INOPHA(IS,I1)
              XB=COD(IB1,1)
              YB=COD(IB1,2)
              ZB=COD(IB1,3)
              IB2=INOPHA(IS,1)
              XC=COD(IB2,1)
              YC=COD(IB2,2)
              ZC=COD(IB2,3)
           END IF
!          �x�N�g��
           X1=-(XA-XB)
           Y1=-(YA-YB)
           Z1=-(ZA-ZB)
           X2=-(XA-XC)
           Y2=-(YA-YC)
           Z2=-(ZA-ZC)
!          �@���x�N�g��
!          NORMAL VECTOR
           A=Y1*Z2-Z1*Y2
           B=Z1*X2-X1*Z2
           C=X1*Y2-Y1*X2
!          HH2���O�łȂ���΁C�ʐς����Q�̃x�N�g���ɂȂ�
!          ���̂Ƃ��̍��W�𗘗p���Ăa�}�g���b�N�X�����
!          IF HH2 IS NOT ZERO, IT BECOMES TWO VECTORS WITH AREA
!          USE THE COORDINATE, MAKE B MATRIX
           HH2=A*A+B*B+C*C
           IF(HH2.GT.0)THEN
              GOTO 300
           END IF
        END DO
300     CONTINUE
!       �v�x�N�g����������p��FLAG
!       �ʍ\���v�f�P�̏d�S�͕��ʂv�x�N�g�����猩�āC�����ɂ���悤�ɂ���
!       FLAG FOR DECISION OF DIRECTION OF W VECTOR
!       SET THE DIRECTION TO MAKE ELEMENT 1 BE UNDER THE W VECTOR DIRECTION

!       AX+BY+CZ+D=0�Ƃ����
!       IF AX+BY+CZ+D=0
        D=-A*XA-B*YA-C*ZA
!       �ʍ\���v�f�P�̏d�S���W
!       CENTER OF GRAVITY OF ELEMENT 1
        IA=IELPHA(IS,1)
        XG1=XYZ(IA,1)
        YG1=XYZ(IA,2)
        ZG1=XYZ(IA,3)
        FLAG=A*XG1+B*YG1+C*ZG1+D
!         FLAG<0�Ȃ�v�͂��̂܂�
!         FLAG>0�Ȃ�v�́@�w-1
!         IF FLAG<0, W DOES NOT CHANGE
!         IF FLAG>0, W*-1
        IF(FLAG.LT.0)THEN
           IFLAG=0
        ELSE 
           IFLAG=-1
        END IF
!       CALCULATION OF R MATRIX
        AUX=XB-XA
        AUY=YB-YA
        AUZ=ZB-ZA
        ACX=XC-XA
        ACY=YC-YA
        ACZ=ZC-ZA
        AWX=AUY*ACZ-AUZ*ACY
        AWY=AUZ*ACX-AUX*ACZ
        AWZ=AUX*ACY-AUY*ACX
!       CANGHE THE DIRECITON OF W DEPENDS ON IFLAG
        IF(IFLAG.EQ.-1)THEN
           AWX=-AWX
           AWY=-AWY
           AWZ=-AWZ
        END IF
        AVX=AWY*AUZ-AWZ*AUY
        AVY=AWZ*AUX-AWX*AUZ
        AVZ=AWX*AUY-AWY*AUX
        XLU=SQRT(AUX**2+AUY**2+AUZ**2)
        XLV=SQRT(AVX**2+AVY**2+AVZ**2)
        XLW=SQRT(AWX**2+AWY**2+AWZ**2)
        EUX=AUX/XLU
        EUY=AUY/XLU
        EUZ=AUZ/XLU
        EVX=AVX/XLV
        EVY=AVY/XLV
        EVZ=AVZ/XLV
        EWX=AWX/XLW
        EWY=AWY/XLW
        EWZ=AWZ/XLW
!       CALCULATION OF B MATRIX
!       CENTER OF COORDINATE OF FACE
        XX=CPHASE(IS,1)
        YY=CPHASE(IS,2)
        ZZ=CPHASE(IS,3)
!       CENTER OF COORDINATE OF ELEMENT 1
        IA=IELPHA(IS,1)
        XG1=XYZ(IA,1)
        YG1=XYZ(IA,2)
        ZG1=XYZ(IA,3)
!       CENTER OF COORDINATE OF ELEMENT 2
        IA=IELPHA(IS,2)
        XG2=XYZ(IA,1)
        YG2=XYZ(IA,2)
        ZG2=XYZ(IA,3)
!       MAKE ARRAY CWORK FOR CALCULATION OF FORCE 
        CWORK(IS,1,1)=EUX
        CWORK(IS,1,2)=EVX
        CWORK(IS,1,3)=EWX
        CWORK(IS,2,1)=EUY
        CWORK(IS,2,2)=EVY
        CWORK(IS,2,3)=EWY
        CWORK(IS,3,1)=EUZ
        CWORK(IS,3,2)=EVZ
        CWORK(IS,3,3)=EWZ
        CWORK(IS,4,1)=XX-XG1
        CWORK(IS,4,2)=YY-YG1
        CWORK(IS,4,3)=ZZ-ZG1
        CWORK(IS,5,1)=XX-XG2
        CWORK(IS,5,2)=YY-YG2
        CWORK(IS,5,3)=ZZ-ZG2
!       CALCULATION OF B MATRIX
        BMAT(IS,1,1)=-EUX
        BMAT(IS,2,1)=-EVX
        BMAT(IS,3,1)=-EWX
        BMAT(IS,1,2)=-EUY
        BMAT(IS,2,2)=-EVY
        BMAT(IS,3,2)=-EWY
        BMAT(IS,1,3)=-EUZ
        BMAT(IS,2,3)=-EVZ
        BMAT(IS,3,3)=-EWZ
        BMAT(IS,1,4)=-(EUY*(ZZ-ZG1)-EUZ*(YY-YG1))
        BMAT(IS,2,4)=-(EVY*(ZZ-ZG1)-EVZ*(YY-YG1))
        BMAT(IS,3,4)=-(EWY*(ZZ-ZG1)-EWZ*(YY-YG1))
        BMAT(IS,1,5)=-(EUZ*(XX-XG1)-EUX*(ZZ-ZG1))
        BMAT(IS,2,5)=-(EVZ*(XX-XG1)-EVX*(ZZ-ZG1))
        BMAT(IS,3,5)=-(EWZ*(XX-XG1)-EWX*(ZZ-ZG1))
        BMAT(IS,1,6)=-(EUX*(YY-YG1)-EUY*(XX-XG1))
        BMAT(IS,2,6)=-(EVX*(YY-YG1)-EVY*(XX-XG1))
        BMAT(IS,3,6)=-(EWX*(YY-YG1)-EWY*(XX-XG1))      

        BMAT(IS,1,7)=EUX
        BMAT(IS,2,7)=EVX
        BMAT(IS,3,7)=EWX
        BMAT(IS,1,8)=EUY
        BMAT(IS,2,8)=EVY
        BMAT(IS,3,8)=EWY
        BMAT(IS,1,9)=EUZ
        BMAT(IS,2,9)=EVZ
        BMAT(IS,3,9)=EWZ
        BMAT(IS,1,10)=(EUY*(ZZ-ZG2)-EUZ*(YY-YG2))
        BMAT(IS,2,10)=(EVY*(ZZ-ZG2)-EVZ*(YY-YG2))
        BMAT(IS,3,10)=(EWY*(ZZ-ZG2)-EWZ*(YY-YG2))
        BMAT(IS,1,11)=(EUZ*(XX-XG2)-EUX*(ZZ-ZG2))
        BMAT(IS,2,11)=(EVZ*(XX-XG2)-EVX*(ZZ-ZG2))
        BMAT(IS,3,11)=(EWZ*(XX-XG2)-EWX*(ZZ-ZG2))
        BMAT(IS,1,12)=(EUX*(YY-YG2)-EUY*(XX-XG2))
        BMAT(IS,2,12)=(EVX*(YY-YG2)-EVY*(XX-XG2))
        BMAT(IS,3,12)=(EWX*(YY-YG2)-EWY*(XX-XG2)) 
     END IF
  END DO

  RETURN
END SUBROUTINE BMATRIX
!--------------------------------------------------------------------------------------------
SUBROUTINE MKARRAY(NNN,IK11,ia,ja,a,XK11,NK11)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::NNN
  INTEGER,INTENT(IN)::IK11(:),NK11(:,:)
  INTEGER,INTENT(OUT)::ia(:),ja(:)
  REAL(KIND=DP),INTENT(IN)::XK11(:,:)
  REAL(KIND=DP),INTENT(OUT)::a(:)

  INTEGER I,II,J,JJ

  II = 0
  IA(1) = 1

  DO I = 1,NNN
     DO J = 1,IK11(I)
        JJ = NK11(I,J)
        IF (JJ >= I)THEN
           II = II+1
           ja(II) = JJ
           a(II) = XK11(I,J)
        ENDIF
     END DO
     ia(I+1) = II+ia(1)
  END DO

  RETURN
END SUBROUTINE MKARRAY
!--------------------------------------------------------------------------------------------
SUBROUTINE FORCE1(NELE,IK12,NK12,XK12,STDEF1,STDEF1BI,B1,NNN,U,ISTEPBI,ISTEP)

  IMPLICIT NONE

  INTEGER,INTENT(IN)::NELE,NNN,ISTEPBI,ISTEP
  INTEGER,INTENT(IN)::IK12(:),NK12(:,:)
  REAL(KIND=DP),INTENT(IN)::XK12(:,:),STDEF1(:),STDEF1BI(:)
  REAL(KIND=DP),INTENT(OUT)::B1(:),U(:)

  INTEGER I,J,AA,BB
  REAL(KIND=DP) XX
  !   SET 0 IN U
  IF(ISTEP.EQ.1)THEN
     U=0.D0
  ENDIF
  !   ITERATION=1
  !   MAKE F OF F=KU
  !   f=-XK12*STDEF
  B1=0.D0
  DO I = 1,NNN
     AA=IK12(I)
     IF(AA.ne.0)THEN
        DO J = 1,AA
           BB=NK12(I,J)
           XX=XK12(I,J)
           IF(ISTEP.EQ.1)THEN
              B1(I)=B1(I)-XX*STDEF1(BB)
           ELSE
              B1(I)=B1(I)-XX*STDEF1BI(BB)
           ENDIF
        END DO
     ENDIF
  END DO

  RETURN
END SUBROUTINE FORCE1
!---------------------------------------------------------------------------
SUBROUTINE ACTUALFORCE1(NELE,NPHASE,KIND0,IELPHA,FA,DDF,U,DU1,IBANE,ISKY,AREA,&
     HH,BANE,STRES,STRA,PRO,BMAT,CWORK,RN,SN,IOUT4,IOUT5,TEN,&
     MTYPE,JYD,ISTEP,JCRA,HIS,HIST,SHIS,SHIST,XYZ,HISP,HISN,HISTP,HISTN)
      
  IMPLICIT NONE
  INTEGER,INTENT(IN)::NELE,NPHASE,KIND0(:),IELPHA(:,:),IBANE(:),ISKY(:),IOUT4,IOUT5,MTYPE,ISTEP
  REAL(KIND=DP),INTENT(IN)::PRO(:,:),BMAT(:,:,:),DU1(:),HH(:,:),BANE(:,:),AREA(:),XYZ(:,:),&
       CWORK(:,:,:),TEN(:)
  INTEGER,INTENT(OUT)::JCRA(:,:)
  REAL(KIND=DP),INTENT(OUT)::FA(:),DDF(:),STRES(:,:),STRA(:,:)
  INTEGER,INTENT(INOUT)::JYD(:,:)
  REAL(KIND=DP),INTENT(INOUT)::U(:),HIS(:,:),HIST(:,:),SHIS(:),SHIST(:),HISP(:,:),HISN(:,:),HISTP(:,:),HISTN(:,:)
  REAL(KIND=DP),INTENT(OUT)::RN,SN
      
  INTEGER I1,IS

  !     CALCULATION OF STRAIN, STRESS AND UNBALANCED FORCE AT ITERATION=1
  !     SET 0 IN UNBALANCED FORCE
  FA = 0.D0
  DDF =0.D0
  !     DISPLACEMENT
    DO I1=1,6*NELE
            U(I1)=U(I1)+1*DU1(I1)     
    END DO
  !     �Ђ��݁C���́C�s�ލ����͂̌v�Z
  !     CALCULATION OF STRAIN, STRESS AND UNBALANCED FORCE
  DO IS=1,NPHASE
     IF(IBANE(IS).EQ.0)THEN
        CALL STRAIN(NELE,STRA,IELPHA,U,BMAT,HH,IS,istep)
        CALL STRESS(STRA,STRES,BANE,PRO,KIND0,IELPHA,IS,TEN,MTYPE,JYD,ISTEP,HH,&
             JCRA,HIS,HIST,SHIS,SHIST,XYZ,HISP,HISN,HISTP,HISTN)
        CALL FORCE(NELE,STRES,FA,CWORK,AREA,IELPHA,IS)
     END IF
  END DO
!     SUMMATION OF SQUARE OF APPLIED FORCE     RN
!     SUMMATION OF SQUARE OF UNBALANCED FORCE  SN
  RN=0
  SN=0
  DO I1=1,6*NELE
     IF(ISKY(I1).GT.0)THEN
        RN=RN+FA(I1)*FA(I1)
     ELSE
        DDF(I1)=FA(I1)
        SN=SN+FA(I1)*FA(I1)
     END IF
  END DO
      
  RETURN
END SUBROUTINE ACTUALFORCE1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE STRAIN(NELE,STRA,IELPHA,U,BMAT,HH,IS,istep)

  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::NELE,IS,istep
  INTEGER,INTENT(IN)::IELPHA(:,:)
  REAL(KIND=DP),INTENT(IN)::U(:),BMAT(:,:,:),HH(:,:)
  REAL(KIND=DP),INTENT(OUT)::STRA(:,:)
  
  INTEGER I1,I2,II,ME,i4
  REAL(KIND=DP)::UU(12)
  REAL(KIND=DP) WW, istrain
  
!     �Ђ��݂̌v�Z
!     CALCULATION OF STRAIN
!     �ψʂ��i�[
!     STORE DISPLACEMENT
      II=0
	  DO 10 I1=1,2
	  ME=IELPHA(IS,I1)
	    DO 20 I2=1,6
		II=II+1
		UU(II)=U((ME-1)*6+I2)
20      CONTINUE
10    CONTINUE

!     �Ђ��݂̌v�Z
!     CALCULATION OF STRAIN
      DO 30 I1=1,3
	  WW=0.0
	    DO 40 I2=1,12
		WW=WW+BMAT(IS,I1,I2)*UU(I2)
40      CONTINUE
      STRA(IS,I1)=WW/(HH(IS,1)+HH(IS,2))
30    CONTINUE
!add on 13112015
!     �Ђ��݂̌v�Z
!     CALCULATION OF NORMAL STRAIN
      WW=0.0
	    DO 50 I2=1,12
		WW=WW+BMAT(IS,3,I2)*UU(I2)
50      CONTINUE
       STRA(IS,3)=WW/(HH(IS,1)+HH(IS,2))        
	  RETURN
END SUBROUTINE STRAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE STRESS(STRA,STRES,BANE,PRO,KIND0,IELPHA,IS, &                     
     TEN,MTYPE,JYD,ISTEP,HH,JCRA,HIS,HIST,SHIS,SHIST,XYZ,HISP,HISN,HISTP,HISTN)

      IMPLICIT NONE
      
      INTEGER,INTENT(IN)::IS,MTYPE,ISTEP
      INTEGER,INTENT(IN)::KIND0(:),IELPHA(:,:)
      REAL(KIND=DP),INTENT(IN)::STRA(:,:),BANE(:,:),PRO(:,:),TEN(:),HH(:,:),XYZ(:,:)
      INTEGER,INTENT(INOUT)::JYD(:,:)
      INTEGER,INTENT(OUT)::JCRA(:,:)
      REAL(KIND=DP),INTENT(INOUT)::STRES(:,:),HIS(:,:),HIST(:,:),SHIS(:),SHIST(:),HISP(:,:),HISN(:,:),HISTP(:,:),HISTN(:,:)
      
      INTEGER ISEI,IST1
      REAL(KIND=DP) FAI0,CC0,FTM,FTMI,XA1,XA2,XA3,XA4,XXA1,XMAXCRACK,XMAXCRACKI,XISEI,&
        SS1,SS2,SN,SS,STS,STN,A,XA,XB,XHH,CRACK,B,FAI,CC,XXA,XXB,XXC,XXD,SSH,DELTA,BBB,BB1,&
        XMAXCRACKII,XMAXCRACKB,restra,restre,SSTA,SSTR,STRAINHD
      
    
      STRES(IS,1)=BANE(IS,1)*STRA(IS,1)
	  STRES(IS,2)=BANE(IS,1)*STRA(IS,2)
      STRES(IS,3)=BANE(IS,2)*STRA(IS,3)
      GOTO 100
        
	  
!     �Ђ���      
!     STRAIN
	  SS1=STRA(IS,1)
	  SS2=STRA(IS,2)
	  SS=SQRT(SS1*SS1+SS2*SS2)
	  SN=STRA(IS,3)
!     �e���ŉ��͌v�Z
!     BEHAVE ELASTICALLY
      STS=BANE(IS,1)*SS
	  STN=BANE(IS,2)*SN

!     �Ђъ��ꂪ�����t���b�O
!

100   CONTINUE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!6666  CONTINUE

      RETURN
END SUBROUTINE STRESS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FORCE(NELE,STRES,FA,CWORK,AREA,IELPHA,IS)
      
  IMPLICIT NONE   
  INTEGER,INTENT(IN)::NELE,IS
  INTEGER,INTENT(IN)::IELPHA(:,:)
      
  REAL(KIND=DP),INTENT(IN)::STRES(:,:),CWORK(:,:,:),AREA(:)
  REAL(KIND=DP),INTENT(OUT)::FA(:)
      
  REAL(KIND=DP)::FG(3)
  REAL(KIND=DP) WW,FG14,FG15,FG16,FG24,FG25,FG26
  INTEGER I1,I2,M1,M2

!     �ߓ_�͂��v�Z����
!     CALCULATE FORCE
!     �Ǐ����W����S�̍��W�n�ɕϊ����āC�ʐς��|����
!     CONVERT FROM LOCAL COORDINATE TO GLOBAL COORDINATE, AND MULTIPLY AREA
	  DO 10 I1=1,3
	  WW=0.0
	    DO 20 I2=1,3
		WW=WW+CWORK(IS,I1,I2)*STRES(IS,I2)
20      CONTINUE
      FG(I1)=WW*AREA(IS)
10    CONTINUE

!     ���[�����g�̌v�Z
!     CALCULATION OF MOMENT
      FG14=CWORK(IS,4,3)*FG(2)-CWORK(IS,4,2)*FG(3)
      FG15=CWORK(IS,4,1)*FG(3)-CWORK(IS,4,3)*FG(1)
	  FG16=CWORK(IS,4,2)*FG(1)-CWORK(IS,4,1)*FG(2)
	  FG24=-CWORK(IS,5,3)*FG(2)+CWORK(IS,5,2)*FG(3)
      FG25=-CWORK(IS,5,1)*FG(3)+CWORK(IS,5,3)*FG(1)
	  FG26=-CWORK(IS,5,2)*FG(1)+CWORK(IS,5,1)*FG(2)

      M1=IELPHA(IS,1)
	  M2=IELPHA(IS,2)
      FA((M1-1)*6+1)=FA((M1-1)*6+1)+FG(1)
	  FA((M1-1)*6+2)=FA((M1-1)*6+2)+FG(2)
      FA((M1-1)*6+3)=FA((M1-1)*6+3)+FG(3)
	  FA((M1-1)*6+4)=FA((M1-1)*6+4)+FG14
	  FA((M1-1)*6+5)=FA((M1-1)*6+5)+FG15
      FA((M1-1)*6+6)=FA((M1-1)*6+6)+FG16

      FA((M2-1)*6+1)=FA((M2-1)*6+1)-FG(1)
	  FA((M2-1)*6+2)=FA((M2-1)*6+2)-FG(2)
      FA((M2-1)*6+3)=FA((M2-1)*6+3)-FG(3)
	  FA((M2-1)*6+4)=FA((M2-1)*6+4)+FG24
	  FA((M2-1)*6+5)=FA((M2-1)*6+5)+FG25
      FA((M2-1)*6+6)=FA((M2-1)*6+6)+FG26

      RETURN
    END SUBROUTINE FORCE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE ACTUALFORCE2(NELE,NPHASE,KIND0,IELPHA,FA,DDF,U,DU,IBANE,&
                        ISKY,AREA,HH,BANE,STRES,STRA,PRO,BMAT,CWORK,RN,SN,&
                        TEN,MTYPE,JYD,ISTEP,JCRA,HIS,HIST,SHIS,SHIST,XYZ,HISP,HISN,HISTP,HISTN)
      
IMPLICIT NONE

INTEGER,INTENT(IN)::NELE,NPHASE,MTYPE,ISTEP
INTEGER,INTENT(IN)::IELPHA(:,:),IBANE(:),ISKY(:),KIND0(:)
REAL(KIND=DP),INTENT(IN)::PRO(:,:),HH(:,:),BANE(:,:),CWORK(:,:,:),BMAT(:,:,:),&
                TEN(:),AREA(:),XYZ(:,:)
REAL(KIND=DP),INTENT(OUT)::STRES(:,:),STRA(:,:),FA(:),DDF(:)
REAL(KIND=DP),INTENT(INOUT)::U(:),DU(:),HIS(:,:),HIST(:,:),SHIS(:),SHIST(:),HISP(:,:),HISN(:,:),HISTP(:,:),HISTN(:,:)
INTEGER,INTENT(INOUT)::JYD(:,:),JCRA(:,:)
REAL(KIND=DP),INTENT(OUT)::RN,SN

INTEGER I1,IS

!     CALCULATION OF STRAIN, STRESS AND UNBALANCED FORCE AT ITERATION>1
!     SET 0 IN UNBALANCED FORCE
      DO 10 I1=1,6*NELE
          FA(I1)=0.0
          DDF(I1)=0.0
10    CONTINUE

!     DISPLACEMENT
      DO 20 I1=1,6*NELE
          U(I1)=U(I1)+DU(I1)
20    CONTINUE

!$omp parallel
!$omp do
!     CALCULATION OF STRAIN, STRESS AND UNBALANCED FORCE
      DO 30 IS=1,NPHASE
	    IF(IBANE(IS).EQ.0)THEN
            CALL STRAIN(NELE,STRA,IELPHA,U,BMAT,HH,IS,istep)
            CALL STRESS(STRA,STRES,BANE,PRO,KIND0,IELPHA,IS,TEN,MTYPE,JYD,ISTEP,HH,JCRA,HIS,HIST,SHIS,SHIST,XYZ,HISP,HISN,HISTP,HISTN)
          !  CALL STRAINB(NELE,STRA,IELPHA,U,BMAT,HH,IS,istep)
            CALL FORCE(NELE,STRES,FA,CWORK,AREA,IELPHA,IS)
        END IF
30    CONTINUE
!$omp end do


!     SUMMATION OF SQUARE OF APPLIED FORCE     RN
!     SUMMATION OF SQUARE OF UNBALANCED FORCE  SN
      RN=0
	  SN=0
!$omp do
DO 40 I1=1,6*NELE
    IF(ISKY(I1).GT.0)THEN
        RN=RN+FA(I1)*FA(I1)
    ELSE
        DDF(I1)=FA(I1)
		SN=SN+FA(I1)*FA(I1)	    
    END IF
40 CONTINUE
!$omp end do

   RETURN
END SUBROUTINE ACTUALFORCE2
!-------------------------------------------------------------------------------------------
SUBROUTINE OUTPUT(NELE,ISTEP,NXDISP,NYDISP,NZDISP,IXDISP,IYDISP,IZDISP,IXDISPI,IYDISPI,IZDISPI,&
                KIND0,IBANE,NNOPHA,INOPHA,NPHELEM,IPHELEM,IELPHA,JCRA,JYD,FA,U,XDISP,YDISP,ZDISP,COD,&
                XYZ,STRA,STRES,IOUT6,IOUT11,IOUT12,IOUT13,IOUT14,IOUT15,IOUT16,IOUT17,IOUT18,IOUT19,IOUT20,&
                IOUT21,IOUT22,IOUT23,IOUT24,IOUT25,IOUT26,IOUT27,IREC,IREC1)

IMPLICIT NONE

INTEGER,INTENT(IN)::NELE,ISTEP,NXDISP,NYDISP,NZDISP,IOUT6,IOUT11,IOUT12,IOUT13,IOUT14,IOUT15,IOUT16,&
                    IOUT17,IOUT18,IOUT19,IOUT20,IOUT21,IOUT22,IOUT23,IOUT24,IOUT25,IOUT26,IOUT27
INTEGER,INTENT(IN)::IXDISP(:),IYDISP(:),IZDISP(:),IXDISPI(:),IYDISPI(:),IZDISPI(:),KIND0(:),IBANE(:),&
                    NNOPHA(:),INOPHA(:,:),NPHELEM(:),IPHELEM(:,:),IELPHA(:,:),JCRA(:,:),JYD(:,:)
INTEGER,INTENT(INOUT)::IREC,IREC1
REAL(KIND=DP),INTENT(IN)::FA(:),U(:),XDISP(:),YDISP(:),ZDISP(:),COD(:,:),XYZ(:,:),STRA(:,:),STRES(:,:)
REAL(KIND=DP) FTO,DF,SDISP,FTOR,SDISPR,FTOL,SDISPL
REAL XX
INTEGER I1,I3,IXX,ITOTAL,IMOR,IAGG,IINT,ITTOTAL,ITMOR,ITAGG,ITINT,ISTOTAL,ISMOR,ISAGG,ISINT,&
        IT2,IT3,IT4,IT5,ITM2,ITM3,ITM4,ITM5,ITA2,ITA3,ITA4,ITA5,ITI2,ITI3,ITI4,ITI5

!     �׏d�ψʋȐ��̕ۑ�
!     REDORD DISPLACEMENT-LOAD CURVE
!     �x�����������׏d�ψʂ��L�^
!     DISP-LOAD OF TOP SIDE IS RECORED
!     �����̍ډה̉׏d�ψʂ��L�^
!     RECORD DISPLACEMENT-LOAD OF SIDE LOADING BOUNDARY

!     0�X�e�b�v�ڂƂ��ă[�������Ă���
!     0 IS WRITTEN FOR STEP 0
      IF(ISTEP==1)THEN
!	  WRITE(IOUT1,8001) 0,0
	  WRITE(IOUT11,8001) 0,0
	  WRITE(IOUT22,8001) 0,0
	  WRITE(IOUT23,8001) 0,0
	  WRITE(IOUT24,8001) 0,0
	  WRITE(IOUT25,8001) 0,0
	  END IF
!eliminate dispload  the !! mean it already ! before if put !
!     ��������
!     LONGITUDINAL DIRECTION
!     �e�X�e�b�v�ł�
!     EACH STEP
!     �׏d
!     LOAD
      FTO=0
	  DO 180 I1=1,NYDISP   !change from Y loading
	  DF=FA((IYDISP(I1)-1)*6+2)   !change from Y loading
      FTO=FTO+DF
180   CONTINUE
!     �ψ�
!     DISPLACEMENT
      SDISP=-U((IYDISP(1)-1)*6+2)	     !change from Y loading
!	  WRITE(IOUT1,8001) SDISP,FTO
      WRITE(IOUT11,8001) SDISP,FTO

!     �ډהE��
!     ���k����
!     RIGTH SIDE BOUNDARY
!     COMPRESSION IS POSITIVE
      FTOR=0
	  DO 190 I1=1,NYDISP   !change from Y loading
	    DF=FA((IYDISP(I1)-1)*6+2)   !change from Y loading
	    FTOR=FTOR+DF
190   CONTINUE
!      SDISPR=0
!      DO 200 I1=1,NXDISP
!        SDISPR=-U((IXDISP(I1)-1)*6+1)
      SDISPR=-U((IYDISP(1)-1)*6+2)   !change from Y loading
200   CONTINUE
	  WRITE(IOUT22,8001) SDISPR,FTOR
      WRITE(IOUT23,8001) SDISPR,FTOR

!     �ډה���
!     ���k����
!     LEFT SIDE BOUNDARY
!     COMPRESSION IS POSITIVE
      FTOL=0
	  DO 210 I1=1,NXDISP
	    DF=-FA((IXDISP(I1)-1)*6+3)
	    FTOL=FTOL+DF
210   CONTINUE
      SDISPL=0
      DO 220 I1=1,NXDISP
        SDISPL=U((IXDISP(I1)-1)*6+3)
220   CONTINUE
	  WRITE(IOUT24,8001) SDISPL,FTOL
      WRITE(IOUT25,8001) SDISPL,FTOL
!end eliminate dispload  the !! mean it already ! before if put !


8001  FORMAT(D15.5,',',D15.5)


!     OUTPUT.INFO�ւ̕ۑ�
!     WRITE IN OUTPUT.INFO
!     �P�X�e�b�v�ڂ̏ꍇ�C���f���̏���ۑ�����
!     IN CASE STEP 1, MODEL INFORMATION IS WRITTEN
      IF(ISTEP.EQ.1)THEN
!     �ߓ_��
!     NUMBER OF NODE
      IXX=NNODE
!     CALL ICONVERT(IXX) 
      IREC=IREC+1
      WRITE(IOUT6) IXX
!     WRITE(IOUT6,REC=IREC) IXX
!     �v�f��
!     NUMBER OF ELEMENT
      IXX=NELE
!     CALL ICONVERT(IXX) 
      IREC=IREC+1
      WRITE(IOUT6) IXX
!     WRITE(IOUT6,REC=IREC) IXX 
!     �ʐ�
!     NUMBER OF FACE
      IXX=NPHASE
!     CALL ICONVERT(IXX) 
      IREC=IREC+1
      WRITE(IOUT6) IXX
!     WRITE(IOUT6,REC=IREC) IXX
	  
!       �ߓ_���W 
!       COORDINATE OF NODE
        DO 10 I1=1,NNODE
	      DO 10 I2=1,3
		  XX=COD(I1,I2)
!         CALL RCONVERT(XX)
          IREC=IREC+1
          WRITE(IOUT6) XX
!          WRITE(IOUT6,REC=IREC) XX
10      CONTINUE
!       �ʍ\���ߓ_�ԍ�
!       NODE NUMBER COMPOSING FACE
        DO 20 I1=1,NPHASE
	    IXX=NNOPHA(I1)
!       CALL ICONVERT(IXX) 
        IREC=IREC+1
        WRITE(IOUT6) IXX
!       WRITE(IOUT6,REC=IREC) IXX
	      DO 20 I2=1,NNOPHA(I1)
		  IXX=INOPHA(I1,I2)
!         CALL ICONVERT(IXX) 
          IREC=IREC+1
          WRITE(IOUT6) IXX
!         WRITE(IOUT6,REC=IREC) IXX
20      CONTINUE	  
!       �v�f�\���ʔԍ�
!       FACE NUMBER COMPOSING ELEMENT
        DO 30 I1=1,NELE
	    IXX=NPHELEM(I1)
!       CALL ICONVERT(IXX) 
        IREC=IREC+1
        WRITE(IOUT6) IXX
!       WRITE(IOUT6,REC=IREC) IXX		 
	      DO 30 I2=1,NPHELEM(I1)
		  IXX=IPHELEM(I1,I2)
!         CALL ICONVERT(IXX) 
          IREC=IREC+1
          WRITE(IOUT6) IXX
!         WRITE(IOUT6,REC=IREC) IXX
30      CONTINUE
!       �ʍ\���v�f�ԍ�
!       ELEMENT NUMBER COMPOSING OF FACE
        DO 40 I1=1,NPHASE
        IXX=IELPHA(I1,1)
!       CALL ICONVERT(IXX) 
        IREC=IREC+1
        WRITE(IOUT6) IXX
!       WRITE(IOUT6,REC=IREC) IXX
	    IXX=IELPHA(I1,2)
!       CALL ICONVERT(IXX) 
        IREC=IREC+1
        WRITE(IOUT6) IXX
!       WRITE(IOUT6,REC=IREC) IXX
40      CONTINUE
!       �v�f��ޔԍ�
!       ELEMENT KIND
        DO 50 I1=1,NELE
	    IXX=KIND0(I1)
!       CALL ICONVERT(IXX) 
        IREC=IREC+1
        WRITE(IOUT6) IXX
!       WRITE(IOUT6,REC=IREC) IXX
50      CONTINUE
!       �o�lFLAG
!       FLAG OF SPRING
        DO 80 I1=1,NPHASE
        IXX=IBANE(I1)
!       CALL ICONVERT(IXX) 
        IREC=IREC+1
        WRITE(IOUT6) IXX
!       WRITE(IOUT6,REC=IREC) IXX
80      CONTINUE
!       �v�f�d�S���W
!       CENTER OF GRAVITY OF ELEMENT
        DO 90 I1=1,NELE
	      DO 90 I2=1,3
		  XX=XYZ(I1,I2)
!         CALL RCONVERT(XX)
          IREC=IREC+1
          WRITE(IOUT6) XX
!      WRITE(IOUT6,REC=IREC) XX
90      CONTINUE
!     �ŏI�X�e�b�v
!     FINAL STEP
      IXX=IFSTEP
!     CALL ICONVERT(IXX) 
      IREC=IREC+1
      WRITE(IOUT6) IXX
!     WRITE(IOUT6,REC=IREC) IXX

      END IF

!     �e�X�e�b�v�ɑ΂��ĕۑ�
!     RECORD AT EACH STEP
!     �׏d�ƕψʂ��L�^
!     REDORD DISPLACEMENT-LOAD CURVE     
!     �㑤�ډה���̈��k�ɂ����Ή����Ă��Ȃ�
!     TOP BOUNDARY IS SELECTED (COMPRESSION IS POSITIVE)
!     �K�v�ȏꍇ�͕ύX
!     MODIFICATION IS NECESSARY FOR OTHER CASES
      XX=SDISP
!     CALL RCONVERT(XX)
      IREC=IREC+1
      WRITE(IOUT6) XX
!     WRITE(IOUT6,REC=IREC) XX
	  XX=FTO
!     CALL RCONVERT(XX)
      IREC=IREC+1
      WRITE(IOUT6) XX
!     WRITE(IOUT6,REC=IREC) XX

!     �e�v�f�̕ψ�
!     DISPLACEMENT OF ELEMENT
      DO 60 I1=1,6*NELE
      XX=U(I1)
!     CALL RCONVERT(XX)
      IREC=IREC+1
      WRITE(IOUT6) XX
!     WRITE(IOUT6,REC=IREC) XX
60    CONTINUE
!     �Ђ��݂Ɖ��͂̕ۑ�
!     RECORD STRAIN AND STRESS
      DO 70 I1=1,NPHASE
	    DO 70 I2=1,3
        XX=STRA(I1,I2)
!       CALL RCONVERT(XX)
        IREC=IREC+1
        WRITE(IOUT6) XX
!       WRITE(IOUT6,REC=IREC) XX
		XX=STRES(I1,I2)
!       CALL RCONVERT(XX)
        IREC=IREC+1
        WRITE(IOUT6) XX
!       WRITE(IOUT6,REC=IREC) XX
70    CONTINUE


!     ���ޕ\�ʂ̂Ђ��݂Ɖ��͂̋L�^
!     STRAIN AND STRESS OF INTERFACE
!     �P�X�e�b�v�ڂ̏ꍇ
!     AT STEP 1
!      �ʂ̍��W���L�^
!     RECORD COORDINATE OF FACE
      I2=0
	  IF(ISTEP.EQ.1)THEN
	    DO 140 I1=1,NPHASE
		  IF(IBANE(I1).EQ.0)THEN
            IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.10.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.11.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.12.OR. &
                KIND0(IELPHA(I1,1)).EQ.10.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.11.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.12.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                    KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN
            I2=I2+1
		    END IF
		  END IF
140     CONTINUE
!     �E�ʖʐ�
!     NUMBER OF INTERFACE
      IXX=I2
      IREC1=IREC1+1
!      WRITE(IOUT17,REC=IREC1) IXX
       WRITE(IOUT17) IXX

        DO 150 I1=1,NPHASE
		  IF(IBANE(I1).EQ.0)THEN
             IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.10.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.11.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.12.OR. &
                KIND0(IELPHA(I1,1)).EQ.10.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.11.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.12.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                    KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN
!           �搔
!           NUMBER OF NODE
            IXX=NNOPHA(I1)
            IREC1=IREC1+1
!      WRITE(IOUT17,REC=IREC1) IXX
      WRITE(IOUT17) IXX
			  DO 160 I2=1,NNOPHA(I1)
			  I3=INOPHA(I1,I2)
!             �ߓ_X���W
!             X COORDINATE
              XX=COD(I3,1)
              IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
!             �ߓ_Y���W
!             Y COORDINATE
              XX=COD(I3,2)
              IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
!             �ߓ_Z���W
!             Z COORDINATE
              XX=COD(I3,3)
              IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
160           CONTINUE
            END IF
          END IF
150     CONTINUE
      END IF
!     ��������S�X�e�b�v�ɂ���
!     FOR ALL STEP
      DO 170 I1=1,NPHASE
	    IF(IBANE(I1).EQ.0)THEN
!         �E�ʂ̏ꍇ�͂Ђ��݂Ɖ��͂��L�^
!         IN CASE INTERFACE, RECORD STRAIN AND STRESS
!         �Ō�ɂЂъ��ꔭ���t���b�O�̋L�^
!         CRACK FLAG IS RECORDED FINALLY
           IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.10.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.11.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.12.OR. &
                KIND0(IELPHA(I1,1)).EQ.10.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.11.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.12.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                    KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN
          XX=STRA(I1,1)
          IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
		  XX=STRA(I1,2)
          IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
		  XX=STRA(I1,3)
          IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
		  XX=STRES(I1,1)
          IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
		  XX=STRES(I1,2)
          IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
		  XX=STRES(I1,3)
          IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) XX
              WRITE(IOUT17) XX
		  
!         ���E�o�l�̏��
!         CONDITION OF SPRING 
!         �j��Ȃ��@�@�@�@0
!         ������Ŕj��@�@1
!         ����f��Ŕj��@6
!         NO DAMAGE       0
!         TENSILE FAILURE 1
!         SHEAR FAILURE   		  
		    IF(JYD(I1,1).EQ.0)THEN
			IXX=0
			ELSE IF(JCRA(I1,5).EQ.1)THEN  !0.03mm
			IXX=5
			ELSE IF(JCRA(I1,4).EQ.1)THEN  !0.01mm
			IXX=4
			ELSE IF(JCRA(I1,3).EQ.1)THEN  !0.005mm
			IXX=3
			ELSE IF(JCRA(I1,2).EQ.1)THEN  !0.002mm
			IXX=2
			ELSE IF(JCRA(I1,1).EQ.1)THEN  !�Ђъ��ꂠ��  WITH CRACK
			IXX=1
			ELSE 
			IXX=6
			END IF
          IREC1=IREC1+1
!             WRITE(IOUT17,REC=IREC1) IXX
              WRITE(IOUT17) IXX
          END IF
        END IF
170   CONTINUE



!     �j�󂵂��o�l�̋L�^
!     DAMAGED SPRING 
      ITOTAL=0
      IMOR=0
	  IAGG=0
	  IINT=0
	  ITTOTAL=0
	  ITMOR=0
	  ITAGG=0
	  ITINT=0
	  ISTOTAL=0
	  ISMOR=0
	  ISAGG=0
	  ISINT=0
	  IT2=0
	  IT3=0
	  IT4=0
	  IT5=0
	  ITM2=0
	  ITM3=0
	  ITM4=0
	  ITM5=0
	  ITA2=0
	  ITA3=0
	  ITA4=0
	  ITA5=0
	  ITI2=0
	  ITI3=0
	  ITI4=0
	  ITI5=0




      IF(ISTEP.EQ.1)THEN
      WRITE(IOUT12,8003)
	  WRITE(IOUT12,8004) 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
      WRITE(IOUT13,8002) 1
	  WRITE(IOUT13,8002) 0
	  WRITE(IOUT14,8002) 1
	  WRITE(IOUT14,8002) 0
	  WRITE(IOUT15,8002) 1
	  WRITE(IOUT15,8002) 0
	  WRITE(IOUT16,8002) 1
	  WRITE(IOUT16,8002) 0
	  WRITE(IOUT18,8002) 1
	  WRITE(IOUT18,8002) 0
	  WRITE(IOUT19,8002) 1
	  WRITE(IOUT19,8002) 0
	  WRITE(IOUT20,8002) 1
	  WRITE(IOUT20,8002) 0
      WRITE(IOUT21,8002) 1
	  WRITE(IOUT21,8002) 0
      WRITE(IOUT26,8002) 1
	  WRITE(IOUT26,8002) 0
      WRITE(IOUT27,8002) 1
	  WRITE(IOUT27,8002) 0


	  ELSE

      WRITE(IOUT18,8002) ISTEP
	  WRITE(IOUT13,8002) ISTEP
	  WRITE(IOUT14,8002) ISTEP
	  WRITE(IOUT15,8002) ISTEP
	  WRITE(IOUT16,8002) ISTEP
	  WRITE(IOUT19,8002) ISTEP
	  WRITE(IOUT20,8002) ISTEP
	  WRITE(IOUT21,8002) ISTEP
	  WRITE(IOUT26,8002) ISTEP
	  WRITE(IOUT27,8002) ISTEP
	   
	    DO 130 I1=1,NPHASE
          IF(JYD(I1,1).EQ.1)THEN
		  ITOTAL=ITOTAL+1
            IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN
            IMOR=IMOR+1
            ELSE IF(KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.&
            KIND0(IELPHA(I1,1)).EQ.10.AND.KIND0(IELPHA(I1,2)).EQ.10.OR.&
            KIND0(IELPHA(I1,1)).EQ.11.AND.KIND0(IELPHA(I1,2)).EQ.11.OR.&
            KIND0(IELPHA(I1,1)).EQ.12.AND.KIND0(IELPHA(I1,2)).EQ.12)THEN
			IAGG=IAGG+1
			ELSE IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.10.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.11.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.12.OR. &
                KIND0(IELPHA(I1,1)).EQ.10.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.11.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.12.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                    KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN
            IINT=IINT+1
			END IF
		  
		    IF(JCRA(I1,1).EQ.1)THEN
			ITTOTAL=ITTOTAL+1
              IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN
              ITMOR=ITMOR+1
			    IF(JCRA(I1,2).EQ.1)THEN
				IT2=IT2+1
				ITM2=ITM2+1
				END IF
				IF(JCRA(I1,3).EQ.1)THEN
				IT3=IT3+1
				ITM3=ITM3+1
				END IF
				IF(JCRA(I1,4).EQ.1)THEN
				IT4=IT4+1
				ITM4=ITM4+1
				END IF
				IF(JCRA(I1,5).EQ.1)THEN
				IT5=IT5+1
				ITM5=ITM5+1
				END IF
					  
			  ELSE IF(KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.2.OR.&
            KIND0(IELPHA(I1,1)).EQ.10.AND.KIND0(IELPHA(I1,2)).EQ.10.OR.&
            KIND0(IELPHA(I1,1)).EQ.11.AND.KIND0(IELPHA(I1,2)).EQ.11.OR.&
            KIND0(IELPHA(I1,1)).EQ.12.AND.KIND0(IELPHA(I1,2)).EQ.12)THEN
			  ITAGG=ITAGG+1
			    IF(JCRA(I1,2).EQ.1)THEN
				IT2=IT2+1
				ITA2=ITA2+1
				END IF
				IF(JCRA(I1,3).EQ.1)THEN
				IT3=IT3+1
				ITA3=ITA3+1
				END IF
				IF(JCRA(I1,4).EQ.1)THEN
				IT4=IT4+1
				ITA4=ITA4+1
				END IF
				IF(JCRA(I1,5).EQ.1)THEN
				IT5=IT5+1
				ITA5=ITA5+1
				END IF
			  			  
			  ELSE IF(KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.2.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.10.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.11.OR. &
                KIND0(IELPHA(I1,1)).EQ.1.AND.KIND0(IELPHA(I1,2)).EQ.12.OR. &
                KIND0(IELPHA(I1,1)).EQ.10.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.11.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                KIND0(IELPHA(I1,1)).EQ.12.AND.KIND0(IELPHA(I1,2)).EQ.1.OR. &
                    KIND0(IELPHA(I1,1)).EQ.2.AND.KIND0(IELPHA(I1,2)).EQ.1)THEN
              ITINT=ITINT+1

                IF(JCRA(I1,2).EQ.1)THEN
				IT2=IT2+1
				ITI2=ITI2+1
				END IF
				IF(JCRA(I1,3).EQ.1)THEN
				IT3=IT3+1
				ITI3=ITI3+1
				END IF
				IF(JCRA(I1,4).EQ.1)THEN
				IT4=IT4+1
				ITI4=ITI4+1
				END IF
				IF(JCRA(I1,5).EQ.1)THEN
				IT5=IT5+1
				ITI5=ITI5+1
				END IF

			  END IF
            END IF
          
!           �Ђъ��ꕝ�ɒB�����ʔԍ��̋L�^
!           FACE NUMBER REACHING CERTAIN CRACK WIDTH
            IF(JCRA(I1,1).EQ.1)THEN
		    WRITE(IOUT18,8002) I1
		    END IF
            IF(JCRA(I1,2).EQ.1)THEN
		    WRITE(IOUT13,8002) I1
		    END IF
		    IF(JCRA(I1,3).EQ.1)THEN
		    WRITE(IOUT14,8002) I1
		    END IF
		    IF(JCRA(I1,4).EQ.1)THEN
		    WRITE(IOUT15,8002) I1
		    END IF
		    IF(JCRA(I1,5).EQ.1)THEN
		    WRITE(IOUT16,8002) I1
		    END IF
			IF(JCRA(I1,6).EQ.1)THEN
		    WRITE(IOUT19,8002) I1
		    END IF
			IF(JCRA(I1,7).EQ.1)THEN
		    WRITE(IOUT20,8002) I1
		    END IF
			IF(JCRA(I1,8).EQ.1)THEN
		    WRITE(IOUT21,8002) I1
		    END IF	
		    IF(JCRA(I1,9).EQ.1)THEN
		    WRITE(IOUT26,8002) I1
		    END IF
			IF(JCRA(I1,10).EQ.1)THEN
		    WRITE(IOUT27,8002) I1
		    END IF	  

		  END IF

130     CONTINUE
      
      ISTOTAL=ITOTAL-ITTOTAL
	  ISMOR=IMOR-ITMOR
	  ISAGG=IAGG-ITAGG
	  ISINT=IINT-ITINT	  
	  
	  WRITE(IOUT12,8004) ISTEP,ITMOR,ISMOR,ITM2,ITM3,ITM4,ITM5,IMOR,ITAGG,ISAGG,ITA2,ITA3,ITA4,ITA5,IAGG,ITINT,ISINT,ITI2,ITI3,ITI4,ITI5,IINT,ITTOTAL,ISTOTAL,IT2,IT3,IT4,IT5,ITOTAL
      
	  WRITE(IOUT18,8002) 0
	  WRITE(IOUT13,8002) 0
	  WRITE(IOUT14,8002) 0
	  WRITE(IOUT15,8002) 0
	  WRITE(IOUT16,8002) 0
	  WRITE(IOUT19,8002) 0
	  WRITE(IOUT20,8002) 0
	  WRITE(IOUT21,8002) 0
	  WRITE(IOUT26,8002) 0
	  WRITE(IOUT27,8002) 0
	  
	  END IF

!       ELSE IF(KIND0(IELPHA(IS,1)).EQ.1.AND.KIND0(IELPHA(IS,2)).EQ.2.OR. &
!              KIND0(IELPHA(IS,1)).EQ.2.AND.KIND0(IELPHA(IS,2)).EQ.1)THEN
8002  FORMAT(I7)
8003  FORMAT('STEP,MOR-TEN,MOR-SHE,MOR-0.002mm,MOR-0.005mm,MOR-0.01mm,MOR-0.03mm,MOR-TOTAL,AGG-TEN,AGG-SHE,AGG-0.002mm,AGG-0.005mm,AGG-0.01mm,AGG-0.03mm,AGG-TOTAL,INT-TEN,INT-SHE,INT-0.002mm,INT-0.005mm,INT-0.01mm,INT-0.03mm,INT-TOTAL,TEN-TOTAL,SHE-TOTAL,0.002-TOTAL,0.005-TOTAL,0.01-TOTAL,0.03-TOTAL,TOTAL')
8004  FORMAT(I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7,',',I7)

RETURN
END SUBROUTINE OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE RELOAD()

!END SUBROUTINE RELOAD
!--------------------------------------------------------------------------------------------
!SUBROUTINE RESIMULATION()

!END SUBROUTINE RESIMULATION
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
