import os
from MMTK import *
from mbpol import mbpolForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest, forceConstantTest, virialTest
from MMTK.Trajectory import Trajectory, TrajectoryOutput, StandardLogOutput, SnapshotGenerator
from MMTK.ForceFields.ForceField import CompoundForceField
from MMTK.Dynamics import TranslationRemover, RotationRemover
from MMTK_PIGSNormalModeIntegrator import PIGSNormalModeIntegrator, PIGSLangevinNormalModeIntegrator
from MMTK_PIGSCartesianIntegrator import PIGSCartesianIntegrator, PIGSLangevinCartesianIntegrator
from MMTK.Trajectory import Trajectory, TrajectoryOutput
from MMTK.NormalModes import VibrationalModes
from MMTK import Features
from sys import argv
from MMTK.Minimization import SteepestDescentMinimizer
from Scientific import N
from numpy import zeros, array, sqrt, asarray, linalg, correlate, fft, real, size, dot, var, arange
from Scientific.Statistics import mean, standardDeviation

###################################################
############ CHOOSE SYSTEM PARAMETERS #############
###################################################

if len(argv) != 4:
	raise("usgage:"+argv[0]+" <traj> <P> <startno> \n") 


#Currently hard-coded, but could make an argument#
strbeta="p006"
beta=float(0.006/Units.k_B)

#################################################
#############  SET UP SYSTEM ####################
#################################################

traj=argv[1]
P=int(argv[2])
if (P % 2 ==0):
        raise()
tau=beta/float(P)
start=int(argv[3])

traj=str(traj)+str(start)+".nc"

temperature=1./(Units.k_B*beta)

bootstrap_dt=.001*Units.ps
bootstrap_friction = 0.01/bootstrap_dt
universe =InfiniteUniverse()

pos1 = Vector(0.0,0.0,0.0)

universe.addObject(Molecule('spcfw-q', position=pos1))
natoms=universe.numberOfAtoms()
for atom in universe.atomList():
        atom.setNumberOfBeads(P)

universe.addObject(Environment.PathIntegrals(temperature, False))


universe.setForceField(mbpolForceField(universe))

universe.environmentObjectList(Environment.PathIntegrals)[0].include_spring_terms = False
universe._changed(True)


####################################################
######### CALC IMAGINARY TIME CORR FUNCTION ########
####################################################

# This "dist" value is the number of beads needed to converge psi_T to psi_0
# (beta_opt/tau)/2

dist=int((0.002/Units.k_B)/(2.0*tau))   
maxP=P/2-dist

c=zeros((3,maxP),float)
c2=zeros((3,maxP),float)
counter=0

label="P-"+str(P)+"-"+str(start)

c0file=open("/warehouse/mdgschmi/MBpolMonomer/corr-sym-"+label,"w")
c1file=open("/warehouse/mdgschmi/MBpolMonomer/corr-bend-"+label,"w")
c2file=open("/warehouse/mdgschmi/MBpolMonomer/corr-asym-"+label,"w")

nmodes=3*natoms-6

r=zeros((natoms,P,3),float)
trajlength= len(Trajectory(universe,traj))
mean0=0.
mean1=0.
mean2=0.
for i in range(trajlength):

        universe.setFromTrajectory(Trajectory(universe,traj),i)

        r[0]=asarray(universe.atomList()[0].beadPositions())  #Hydrogen 1
        r[1]=asarray(universe.atomList()[1].beadPositions())  #Hydrogen 2
        r[2]=asarray(universe.atomList()[2].beadPositions())  #Oxygen

        bond_down=zeros((3,maxP),float)
        bond_up=zeros((3,maxP),float)
        
        for p in range(maxP):
                #down refers to correlating the function from the middle to the "left" of the path
                # up  refers to correlating the function from the middle to the "right" of the path

                #[0] corresponds to correlating OH1 + OH2 distance : Symmetric Stretch
                bond_down[0][p]  = (linalg.norm(r[2][(P-1)/2-p]-r[0][(P-1)/2-p])+linalg.norm(r[2][(P-1)/2-p]-r[1][(P-1)/2-p]))/2.0
                bond_up[0][p]    = (linalg.norm(r[2][(P-1)/2+p]-r[0][(P-1)/2+p])+linalg.norm(r[2][(P-1)/2+p]-r[1][(P-1)/2+p]))/2.0

                #[1] corresponds to correlating cos(theta) : Bend motion
                bond_down[1][p]  = dot((r[2][(P-1)/2-p]-r[0][(P-1)/2-p]),(r[2][(P-1)/2-p]-r[1][(P-1)/2-p]))/(linalg.norm(r[2][(P-1)/2-p]-r[0][(P-1)/2-p])*linalg.norm(r[2][(P-1)/2-p]-r[1][(P-1)/2-p]))
                bond_up[1][p]    = dot((r[2][(P-1)/2+p]-r[0][(P-1)/2+p]),(r[2][(P-1)/2+p]-r[1][(P-1)/2+p]))/(linalg.norm(r[2][(P-1)/2+p]-r[0][(P-1)/2+p])*linalg.norm(r[2][(P-1)/2+p]-r[1][(P-1)/2+p]))

                #[2] corresponds to correlating OH1-OH2 distance: Anti-Symmetric Stretch
                bond_down[2][p]  = (linalg.norm(r[2][(P-1)/2-p]-r[0][(P-1)/2-p])-linalg.norm(r[2][(P-1)/2-p]-r[1][(P-1)/2-p]))/2.0
                bond_up[2][p]    = (linalg.norm(r[2][(P-1)/2+p]-r[0][(P-1)/2+p])-linalg.norm(r[2][(P-1)/2+p]-r[1][(P-1)/2+p]))/2.0

        mean0+=mean(bond_down[0])+mean(bond_up[0])
        mean1+=mean(bond_down[1])+mean(bond_up[1])
        mean2+=mean(bond_down[2])+mean(bond_up[2])
        

mean0/=2.0*trajlength
mean1/=2.0*trajlength
mean2/=2.0*trajlength

for i in range(trajlength):

        universe.setFromTrajectory(Trajectory(universe,traj),i)

        r[0]=asarray(universe.atomList()[0].beadPositions())  #Hydrogen 1
        r[1]=asarray(universe.atomList()[1].beadPositions())  #Hydrogen 2
        r[2]=asarray(universe.atomList()[2].beadPositions())  #Oxygen

        bond_down=zeros((3,maxP),float)
        bond_up=zeros((3,maxP),float)
        
        for p in range(maxP):
                #down refers to correlating the function from the middle to the "left" of the path
                # up  refers to correlating the function from the middle to the "right" of the path

                #[0] corresponds to correlating OH1 + OH2 distance : Symmetric Stretch
                bond_down[0][p]  = (linalg.norm(r[2][(P-1)/2-p]-r[0][(P-1)/2-p])+linalg.norm(r[2][(P-1)/2-p]-r[1][(P-1)/2-p]))/2.0
                bond_up[0][p]    = (linalg.norm(r[2][(P-1)/2+p]-r[0][(P-1)/2+p])+linalg.norm(r[2][(P-1)/2+p]-r[1][(P-1)/2+p]))/2.0

                #[1] corresponds to correlating cos(theta) : Bend motion
                bond_down[1][p]  = dot((r[2][(P-1)/2-p]-r[0][(P-1)/2-p]),(r[2][(P-1)/2-p]-r[1][(P-1)/2-p]))/(linalg.norm(r[2][(P-1)/2-p]-r[0][(P-1)/2-p])*linalg.norm(r[2][(P-1)/2-p]-r[1][(P-1)/2-p]))
                bond_up[1][p]    = dot((r[2][(P-1)/2+p]-r[0][(P-1)/2+p]),(r[2][(P-1)/2+p]-r[1][(P-1)/2+p]))/(linalg.norm(r[2][(P-1)/2+p]-r[0][(P-1)/2+p])*linalg.norm(r[2][(P-1)/2+p]-r[1][(P-1)/2+p]))

                #[2] corresponds to correlating OH1-OH2 distance: Anti-Symmetric Stretch
                bond_down[2][p]  = (linalg.norm(r[2][(P-1)/2-p]-r[0][(P-1)/2-p])-linalg.norm(r[2][(P-1)/2-p]-r[1][(P-1)/2-p]))/2.0
                bond_up[2][p]    = (linalg.norm(r[2][(P-1)/2+p]-r[0][(P-1)/2+p])-linalg.norm(r[2][(P-1)/2+p]-r[1][(P-1)/2+p]))/2.0


        bond_down[0]=bond_down[0]-mean0
        bond_down[1]=bond_down[1]-mean1
        bond_down[2]=bond_down[2]-mean2
        
        bond_up[0]=bond_up[0]-mean0
        bond_up[1]=bond_up[1]-mean1
        bond_up[2]=bond_up[2]-mean2

        c0down=zeros(maxP,float)
        c0up=zeros(maxP,float)

        c1down=zeros(maxP,float)
        c1up=zeros(maxP,float)
        
        c2down=zeros(maxP,float)
        c2up=zeros(maxP,float)
        
        for p in range(maxP):
                c0down[p]=bond_down[0][0]*bond_down[0][p]
                c0up[p]=bond_up[0][0]*bond_up[0][p]

                c1down[p]=bond_down[1][0]*bond_down[1][p]
                c1up[p]=bond_up[1][0]*bond_up[1][p]

                c2down[p]=bond_down[2][0]*bond_down[2][p]
                c2up[p]=bond_up[2][0]*bond_up[2][p]
        
        c[0]+=0.5*(c0down+c0up)
        c[1]+=0.5*(c1down+c1up)
        c[2]+=0.5*(c2down+c2up)

        c2[0]+=0.5*(array(c0down)*array(c0down)+array(c0up)*array(c0up))
        c2[1]+=0.5*(array(c1down)*array(c1down)+array(c1up)*array(c1up))
        c2[2]+=0.5*(array(c2down)*array(c2down)+array(c2up)*array(c2up))

c=c/trajlength
c2=c2/trajlength

v0=(c2[0]-(array(c[0]))*(array(c[0])))
v1=(c2[1]-(array(c[1]))*(array(c[1])))
v2=(c2[2]-(array(c[2]))*(array(c[2])))

for i in range(maxP):
        c0file.write(str(i*tau*Units.k_B)+" "+str(c[0][i])+" "+str(v0[i])+"\n")
        c1file.write(str(i*tau*Units.k_B)+" "+str(c[1][i])+" "+str(v1[i])+"\n")
        c2file.write(str(i*tau*Units.k_B)+" "+str(c[2][i])+" "+str(v2[i])+"\n")
        
        
c0file.close()
c1file.close()
c2file.close()
