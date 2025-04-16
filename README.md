# 3D-RBSM
3D Rigid Body Spring Model for Simulating the Mechanical Properties of Concrete
# **Part I: Setting up Fortran Development Environment on Windows Using Visual Studio + Intel oneAPI base Toolkit+ Intel oneAPI HPC Toolkit**

# **1\. Prerequisites**

## **1.1. Visual Studio**

1. Download from: <https://visualstudio.microsoft.com/>
2. Recommended version: **Visual Studio 2019 (Community Edition is sufficient)**
3. During installation, make sure to select:

**Desktop development with C++**

Optionally: **.NET desktop development** (for GUI tools)

## **1.2. Intel oneAPI Base Toolkit**

1. Download from: <https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html>

## **1.3. Intel oneAPI HPC Toolkit**

1. Download from: <https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html>
2. Includes Intel Fortran Compiler (ifx and ifort) and MPI libraries.

Note: After installing Visual Studio, install the Base Toolkit **first**, then the HPC Toolkit.

# **2\. Verifying Installation**

After installation, you should have:

1. Intel oneAPI Command Prompt shortcuts (in Start Menu)

x64 Native Tools Command Prompt for VS 2019

Intel oneAPI Command Prompt for Visual Studio 2019

1. Compilers: ifx.exe, ifort.exe (can be used in both command line and Visual Studio)

# **3\. Compile and Run Fortran Code Using Visual Studio**

## **3.1. Open Visual Studio.**

## **3.2. Create a new Fortran project**

- 1. Go to: File > New > Project
  2. Choose: **Intel Fortran Console Application**
  3. Name your project and solution

## **3.3. Add your .f90 source file (**3D-RBSM-MSM046.f90**)**

- 1. Right-click on **Source Files > Add > Existing Item**
  2. Select your .f90 file.

## **3.4. Find the path of “mkl_pardiso.f90”** 
**For example: C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2024.1\\include**

- 1. Add “mkl_pardiso.f90” to source file (same procedure as adding source file “3D-RBSM-MSM046.f90”)
  2. Program (right-click) > property > general > use complier > IFORT Intel® Fortran Compiler Classic
  3. Program (right-click) > property > Fortran > Additional Include Directories > add the path of “mkl_pardiso.f90” (C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2024.1\\include)

## **3.5. Find the path of “mkl_rt.lib” and “libiomp5md.lib”**  
**For example: “C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2024.1\\lib” and “C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\2024.1\\lib”**

- 1. Program (right-click) > property > Fortran > Linker > General > Additional Include Directories > add the path of “mkl_rt.lib” and “libiomp5md.lib” (C:\\Program Files (x86)\\Intel\\oneAPI\\mkl\\2024.1\\lib C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\2024.1\\lib)
  2. Program (right-click) > property > Fortran > Linker > Input > Additional Dependencies > add “mkl_rt.lib libiomp5md.lib”

## **3.6. Build**

- 1. Build > Build Solution (Ctrl+Shift+B)

## **3.7. Run**

- 1. Debug > Start Without Debugging (Ctrl+F5)

# **Part II: Explanation of the Fortran Program**

**Structure of FORTRAN Program**
![Structure of FORTRAN Program](assets/structure%20of%20RBSM.jpg)

# **1\. READ DATA**

## **1.1. Read model geometrical information from “INDATA1.INFO”**

A model with the size of 100 x 100 x 100 mm is made as the example case (Fig. 1).
![Model size](assets/model%20size.jpg)

Fig. 1 Model size

### **1.1.1. Total number of nodes**

### **1.1.2. Total number of elements**

### **1.1.3. Total number of phases**

The format of the input:
| NUMBER OF NODE |=5958 |
| --- | --- |
| NUMBER OF ELEMENT |=1200 |
| NUMBER OF PHASE |=6955 |

### **1.1.4. Coordinate of node (from the first to the last node): x (change line) y (change line) z (change line)**

The format of the input (the first six nodes for example):

| COORDINATE OF NODE |     |
| --- | --- |
| 43.483223 | \*The x, y, z coordinate of the first node of the model |
| 49.487068 |
| 66.778252 |
| 40.415455 | \*The x, y, z coordinate of the second node of the model |
| 50.288898 |
| 62.514286 |
| 40.923862 | \*The x, y, z coordinate of the third node of the model |
| 50.492710 |
| 61.382435 |
| 6.657995 | \*The x, y, z coordinate of the fourth node of the model |
| 90.561119 |
| 19.351658 |
| 7.162793 | \*The x, y, z coordinate of the fifth node of the model |
| 90.996971 |
| 20.000000 |
| 6.571577 | \*The x, y, z coordinate of the sixth node of the model |
| 71.597572 |
| 23.024349 |

### **1.1.5. Node number composing face**

For one face of each rigid elements, it composes several nodes. In this part, writing the total nodes, and node number of each node composing the face.

The format of the input (the first six faces for example):

| NODE NUMBER COMPOSING FACE |     |
| --- | --- |
| 5_113_370_1959_1960_114 | \*Face 1. The first number is the total node number, after that, the node number of each node composing the face is written one by one. |
| 6_113_371_372_116_115_114 | \*Face 2 |
| 3_113_371_370 | \*Face 3 |
| 8_114_115_363_1815_1816_905_1964_1960 | \*Face 4 |
| 5_115_363_364_117_116 | \*Face 5 |
| 6_116_372_913_912_366_117 | \*Face 6 |

### **1.1.6. Face number composing element**

For each rigid elements, it composes several faces. In this part, writing the total faces, and face number of each face composing the element.

The format of the input (the first six elements for example):

| FACE NUMBER COMPOSING ELEMENT |     |
| --- | --- |
| 18_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18 | \*Element 1. The first number is the total face number, after that, the face number of each face composing the element is written one by one. |
| 17_8_19_20_21_22_23_24_25_26_27_28_29_30_31_32_33_34 | \*Element 2 |
| 12_24_35_36_37_38_39_40_41_42_43_44_45 | \*Element 3 |
| 19_42_46_47_48_49_50_51_52_53_54_55_56_57_58_59_60_61_62_63 | \*Element 4 |
| 15_47_64_65_66_67_68_69_70_71_72_73_74_75_76_77 | \*Element 5 |
| 16_65_78_79_80_81_82_83_84_85_86_87_88_89_90_91_92 | \*Element 6 |

### **1.1.7. Element number composing face**

In the model, several elements may share on face. Therefore, one face connects several elements. In this part, writing the element number connected by the face, and face number.

The format of the input (the first eleven faces for example):

| ELEMENT NUMBER COMPOSING FACE |     |
| --- | --- |
| 1   | \*One element is connected by the face |
| 523 | \*The face number |
| 1   | \*One element is connected by the face |
| 713 | \*The face number |
| 1   | \*One element is connected by the face |
| 513 | \*The face number |
| 1   | \*One element is connected by the face |
| 873 | \*The face number |
| 1   | \*One element is connected by the face |
| 881 | \*The face number |
| 1   | \*One element is connected by the face |
| 723 | \*The face number |
| 1   | \*One element is connected by the face |
| 889 | \*The face number |
| 1   | \*One element is connected by the face |
| 2   | \*The face number |
| 1   | \*One element is connected by the face |
| 714 | \*The face number |
| 1   | \*One element is connected by the face |
| 9   | \*The face number |
| 2   | \*Two elements are connected by the face |
| 881 | \*The face number |

### **1.1.8. Element kind**

Different element is assigned with different kind number (for example, kind number of mortar element is 1, and kind number of aggregate element is 6). The element kind number of the first to the last element is written one by one.

The format of the input (the first six elements):

| ELEMENT KIND NUMBER |     |
| --- | --- |
| 1   | \*Mortar element (the first element) |
| 1   | \*Mortar element (the second element) |
| 1   | \*Mortar element (the third element) |
| 1   | \*Mortar element (the fourth element) |
| 1   | \*Mortar element (the fifth element) |
| 1   | \*Mortar element (the sixth element) |

### **1.1.9. Fixed element condition**

In the RBSM model, displacement is assigned to the elements on the model surface to achieve different boundary conditions. The elements on the model surface may be assigned a displacement or be fixed. Each element has three translational degrees of freedom (x, y, z) and one rotational degree of freedom (θ). First, the total number of elements that have been fixed in the x, y, z directions and for θ should be listed one by one. Then, the fixed element numbers should be listed one by one for the x, y, z directions and θ.

The format of the input:

| NUMBER OF FIXED ELEMENT ON EACH DIRECTION |     |
| --- | --- |
| 200 | \*Total number of elements which have been fixed in x direction |
| 100 | \*Total number of elements which have been fixed in y direction |
| 200 | \*Total number of elements which have been fixed in z direction |
| 200 | \*Total number of elements which have been fixed in |
| FIXED ELEMENT NUMBER ON EACH DIRECTION |     |
| 1001 | \*Each element number for those which have been fixed in x direction; 200 in total. |
| 1002 |
| 1003 |
| …   |
| 1101 | \*Each element number for those which have been fixed in y direction; 100 in total. |
| 1102 |
| 1103 |
| …   |
| 1001 | \*Each element number for those which have been fixed in z direction; 200 in total. |
| 1002 |
| 1003 |
| …   |
| 1001 | \*Each element number for those which have been fixed in (rotation); 200 in total. |
| 1002 |
| 1003 |
| …   |

### **1.1.10. Force displacement condition**

The elements on the model surface may be assigned either a displacement or be fixed. Each element has three translational degrees of freedom (x, y, z) and one rotational degree of freedom (θ). First, the total number of elements that have been given a displacement in the x, y, and z directions should be listed, respectively. Then, for each direction (x, y, and z), the element numbers that have been assigned a displacement should be listed one by one, along with the direction of displacement and the corresponding value.

In this case, the top surface is composed of 100 elements. These elements are given a displacement in the y-direction, while the x, z directions and θ are fixed.

The format of the input:

| NUMBER OF FORCE DISPLACEMENT ELEMENT ON EACH DIRECTION |     |
| --- | --- |
| 0   | \*Total number of elements which have been given displacement in x direction |
| 100 | \*Total number of elements which have been given displacement in y direction |
| 0   | \*Total number of elements which have been given displacement in z direction |
| ELEMENT NUMBER OF FORCE DISPLACEMENT AND THE VALUE |     |
| 1001 | \*Element number |
| 1   | \*Direction of the displacement (1 is in y direction) |
| \-0.01 | \*Displacement value at each calculation step |
| 1002 | \*Element number |
| 1   | \*Direction of the displacement |
| \-0.01 | \*Displacement value at each calculation step |
| 1003 | \*Element number |
| 1   | \*Direction of the displacement |
| \-0.01 | \*Displacement value at each calculation step |
| …   | \*100 groups in total |

### **1.1.11. Final step**

Record the number of calculation steps.

The format of the input:

| FINAL STEP |     |
| --- | --- |
| 200 | \*Simulation will stop at step 200 |

## **1.2. Read material properties from “INDATA2.INFO”**

An RBSM model may contain different types of elements, which are assigned different element kind numbers and have different material properties. For example, in a model that includes aggregate elements and mortar elements, the element kind number for mortar elements is 1, and the element kind number for aggregate elements is 6.

In this case, the model contains only mortar elements; therefore, the element kind number '1' is assigned to the mortar elements. The modulus of elasticity (PRO(1,1)), Poisson's ratio (PRO(1,2)), and tensile strength (PRO(1,3)) are read by the FORTRAN program.

\*Note: In the FORTRAN program, PRO(X,Y) represents the material properties of elements. 'X' is the element kind number, and 'Y' ranges from 1 to 3, representing different properties, including the modulus of elasticity, Poisson's ratio, and tensile strength. The range of 'Y' can be extended according to different conditions.

The format of the input:

| Material properties |     |
| --- | --- |
| Mortar |     |
| Modulus of Elasticity |     |
| 20689.000000 | \*Elastic modulus of mortar elements |
| Poisson Ratio |     |
| 0.180000 | \*Poisson ratio of mortar elements |
| Tensile Strength |     |
| 2.000000 | \*Tensile strength of mortar elements |

**During the simulation, the files "INDATA1.INFO" and "INDATA2.INFO" should be placed in the same folder as the Fortran source file "3D-RBSM-MSM046.f90".**

# **2\. Subroutine INDATA**

## **2.1. Subroutine PHASECENTER**

Calculating the coordinate of the computational point for each face.

Three individual springs including one normal and two shear springs are set at the computational point on the face between two elements. The computational point (_x<sub>cf</sub>_, _y<sub>cf</sub>_, _z<sub>cf</sub>_) is defined (Eq. 1 to Eq. 3),

![eq](assets/eq1%203.jpg)

where, _n_ is the node number composing one face. _x<sub>1</sub>_, _x<sub>2</sub>_, …, _x<sub>n</sub>_ is the _x_ coordinate of each node. _y<sub>1</sub>_, _y<sub>2</sub>_, …, _y<sub>n</sub>_ is the _y_ coordinate of each node. _z<sub>1</sub>_, _z<sub>2</sub>_, …, _z<sub>n</sub>_ is the _z_ coordinate of each node.

## **2.2. Subroutine ELEMCENTER**

Calculating the coordinate of the computational point for each element.

Each element has three transitional and three rotational degrees of freedoms at some point within the element. The computational point (_x<sub>ce</sub>_, _y_<sub>ce</sub>, _z<sub>ce</sub>_) is defined as (Eq. 4 to Eq. 6),
![eq](assets/eq%204%206.jpg)

where, _m_ is the face number composing one element. _x<sub>1</sub>_, _x<sub>2</sub>_, …, _x<sub>m</sub>_ is the _x_ coordinate of the computational point for each face. _y<sub>1</sub>_, _y<sub>2</sub>_, …, _y<sub>m</sub>_ is the _y_ coordinate of the computational point for each face. _z<sub>1</sub>_, _z<sub>2</sub>_, …, _z<sub>m</sub>_ is the _z_ coordinate of the computational point for each face.

## **2.3. Subroutine PHASEAREA**

Calculating the area of each face.

Each face is divided into several triangles from the computational point. The total area of the face is the sum of the areas of each triangle.

## **2.4. Subroutine PERPENDICULAR**

Calculating the normal vector of each face.

## **2.5. Subroutine TMAT**

Calculating the size of displacement matrix K11.

# **3\. Subroutine SPRINGSTIF**

In this spring, the local stiffness of each spring is calculated.

Each face has one normal and two shear springs. _k_<sub>n</sub>, _k_<sub>s1</sub> and _k_<sub>s2</sub> are the normal and shear spring stiffness.

The local stiffness _k<sub>n</sub>,_ _k<sub>s1</sub>_ and _k<sub>s2</sub>_ can be written as (Eq. 7 to Eq. 9),
![eq](assets/eq%207%209.jpg)

where _h<sub>1</sub>_ and _h<sub>2</sub>_ are the length of perpendicular lines from the element computational point to the face where springs are set. _A_ is the area of the face.

_k<sub>nsp</sub>_ and _k<sub>ssp</sub>_ can be calculated under plain strain condition, those are (Eq. 10 and Eq. 11),
![eq](assets/eq%2010%2011.jpg)

_E<sub>elem</sub>_ and _ν<sub>elem</sub>_ are the modulus of elasticity and Poisson’s ratio set on boundary, respectively, which are given by a weighted average of the material properties in two elements according to their perpendicular lengths. Those are (Eq. 12 and Eq. 13),

![eq](assets/eq%2012%2013.jpg)

where _h<sub>1</sub>_ and _h<sub>2</sub>_ are the length of perpendicular lines from the computational point of element1 and 2 to their common face.

The mesoscopic Poisson’s ratio _v<sub>elem</sub>_ is calculated from the macroscopic Poisson’s ratio _v_. The mesoscopic elastic modulus for each element _E<sub>elem</sub>_ is calculated based on the macroscopic elastic modulus _E_ and mesoscopic Poisson’s ratio _v<sub>elem</sub>_ according to following equations (Eq. 14 and Eq. 15),
![eq](assets/eq%2014%2015.jpg)

# **4\. Subroutine SKY**

Recording the boundary conditions for each element.

# **5\. Subroutine KLOCATION**

Recording the non-zero cell for the displacement matrix K11.

# **6\. Subroutine MAKEMATRIX**

Making stiffness matrix for local coordinate system.

## **6.1. Subroutine ELEMATRIX**

### **6.1.1. Subroutine BMATRIX0**

![RBSM elements](assets/RBSM%20elements.jpg)

Fig. 2 RBSM elements

Making the transformation matrix _B_ (from _x_\-_y-z_ to _s1_\-_s2-n_).

The transformation matrix _B_ which is written as,

![B matrix](assets/B%20matrix.jpg)

where _e<sub>ij</sub>_ is direction cosine in _i_ axis on _j_ axis.

### **6.1.2. Subroutine MKELEMAT**

Making local stiffness matrix, and transforming to global coordinate (_k<sub>e</sub>_) using matrix B (Eq. 16).
![eq](assets/eq%2016.jpg)

where,

![eq](assets/eq%2016D.jpg)

in which _k<sub>n</sub>, k<sub>s1</sub>_ and _k<sub>s2</sub>_ are the normal and shear spring stiffness as shown in Eq. 7 to Eq. 9.

### **6.2. Subroutine MKTMAT**

Install stiffness matrix under global coordinate.

# **7\. Subroutine SIMULATION**

## **7.1. Subroutine BMATRIX**

Making transformation matrix _B_ (local coordinate system to global coordinate system, see 7.1.1).

## **7.2. Subroutine MKARRY**

Calculating related parameters for PARDISO matrix solver.

## **7.3. Subroutine FORCE1**

By applying the principle of virtual work, the local equilibrium relation expressed in global coordinate is expressed as,
![eq](assets/eq%2017.jpg)

For each simulation step, under global coordinate is calculated according to Eq. 17 with the given displacement of surface elements.

## **7.4. SMSOLVE1**

After obtain , of each element are calculated by applying the principle of virtual work, as expressed in Eq. 18 for each simulation step,

![eq](assets/eq%2018.jpg)

## **7.5. Subroutine ACTURALFORCE1**

### **7.5.1. Subroutine STRAIN**

At each simulation step, for surface elements which have been assigned displacement, the assigned displacement is directly adopted as the final displacement under global coordinate, while for the elements which have not been given any external displacement, the calculated displacement (Eq. 18) is adopted as the final displacement under global coordinate. The obtained final displacement of each spring under globe coordinate is then converted to local coordinate system (elongation of the springs) using transformation matrix _B_. Then, the strain for springs at each face is calculated based on the elongation value and the original length of each spring.

### **7.5.2. Subroutine STRESS**

The shear and normal stress are calculated.

The behavior of shear and normal springs are simply set as elastic in the program (the stress-strain relationship should be modified according to different conditions and spring types).

The stiffness and are written as shown in Eq. 10 and Eq. 11. The shear stress and normal stress can be expressed as shown in Eq. 19 to Eq. 21.

![eq](assets/eq%2019%2021.jpg)

Where _STRE<sub>1</sub>, STRE <sub>2</sub>_ and _ STRE <sub>3</sub>_ is the stress of the first shear spring, second shear spring, and the normal spring. , and  _STRA<sub>1</sub>, STRA <sub>2</sub>_ and _ STRA <sub>3</sub>_ is the strain of the first shear spring, second shear spring, and the normal spring, respectively.

### **7.5.3. Subroutine FORCE**

Converting stress from local coordinate to global coordinate and multiply area.

## **7.6. Subroutine OUTPUT**

### **7.6.1. Load-displacement information**

The displacement set for each simulation step, as well as the total force of the elements which have been assigned with displacement is recorded in text file “DISP-LOAD(TEXT)” when the displacement is set in _y_ direction. It can be modified according to the direction in which the displacement is set.

### **7.6.2. Stress and displacement information (location) for each element within the model**

All the information for each simulation step was recorded in “OUTPUT.INFO”.

1. The model geometric information is recorded at the first simulation step;
2. The load and displacement for each simulation step;
3. The displacement for each element at each simulation step;
4. The stress and strain for normal and shear springs for each face at each simulation step.

# **8\. Simulation results**

## **8.1. Stress-strain relationship**

The stress-strain curve of the concrete cube under compressive loading is presented in Fig. 3.
![ Stress-strain relationship](assets/stress%20strain.jpg)

Fig. 3 Stress-strain relationship

## **8.2. Stress development of center section**

The stress development of center cross-section of the concrete cube under compressive loading is presented in Fig. 4.
![ Stress-strain relationship](assets/stress.gif)

Stress range: 3 MPa (tension) to 100 MPa (compression)

Fig. 4 Stress-strain relationship
# **References**
1. Nagai, K. et al. (2004). *Mesoscopic Simulation of Failure of Mortar and Concrete by 2D RBSM*. [Journal of Advanced Concrete Technology Vol. 2, No. 3, 359-374, October 2004](https://doi.org/10.3151/jact.2.359) 
2. Nagai, K. et al. (2005). *Mesoscopic simulation of failure of mortar and concrete by 3D RBSM*. [Journal of Advanced Concrete Technology Vol. 3, No. 3, 385-402, October 2005](https://doi.org/10.3151/jact.3.385) 
