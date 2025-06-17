# -*- coding: utf-8 -*-
# 240612

import time

from openmm.app import *
from openmm import *
from openmm import unit
from sys import stdout

ts = time.perf_counter()

# 系の準備
pdb = PDBFile("data/3poz_processed.pdb")  # 前処理済みファイル
forcefield = ForceField("amber14-all.xml","amber14/tip3pfb.xml") # 力場

# System構築に使うパラメータ設定
nonbondedMethod = PME
nonbondedCutoff = 1.0 * unit.nanometers
ewaldErrorTolerance = 5e-4
constraints = HBonds
rigidWater = True
constraintTolerance = 1e-6

# 積分計算に使うパラメータ設定
dt = 0.002 * unit.picoseconds
temperature = 300 * unit.kelvin
friction = 1.0 / unit.picosecond
pressure = 1.0 * unit.atmospheres
barostatInterval = 25

# シミュレーション設定
steps = 500000
equilibrationSteps = 50000

# platform設定
platform = Platform.getPlatformByName("CUDA")
platformProperties = {"Precision":"single"}

# reporter定義
dcdReporter = DCDReporter("result/trajectory.dcd",10000)
dataReporter = StateDataReporter("result/log.txt",1000,totalSteps=steps,step=True,speed=True,progress=True,potentialEnergy=True,temperature=True,separator="\t")
checkpointReporter = CheckpointReporter("result/checkpoint.chk",10000)

# シミュレーション準備
system = forcefield.createSystem(pdb.topology,
                                 nonbondedMethod=nonbondedMethod,
                                 nonbondedCutoff=nonbondedCutoff,
                                 constraints=constraints,
                                 rigidWater=rigidWater,
                                 ewaldErrorTolerance=ewaldErrorTolerance) # PME: 粒子・メッシュ・エバルト法
system.addForce(MonteCarloBarostat(pressure,temperature,barostatInterval)) # 圧力を設定することでNPTアンサンブルにする

# 積分計算設定
integrator = LangevinMiddleIntegrator(temperature,friction,dt)  # ランジュバン動力学に基づいた積分計算
integrator.setConstraintTolerance(constraintTolerance)

# シミュレーション構築
simulation = Simulation(pdb.topology,system,integrator,platform,platformProperties)
simulation.context.setPositions(pdb.positions)

# xml出力設定
with open("result/system.xml",mode="w") as f:
    f.write(XmlSerializer.serialize(system))
with open("result/integrator.xml",mode="w") as f:
    f.write(XmlSerializer.serialize(integrator))

# エネルギー最小化
print("energy minimization start")
simulation.minimizeEnergy()

# 平衡化 → エネルギー最小化後の粒子に初速を与える
print("equilibration start")
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

# reporterの追加
simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0

# シミュレーション
print("simulation start")
simulation.step(steps)

# PDBx/mmCIFを出力する
simulation.saveState("result/final_state.xml")
state = simulation.context.getState(getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions())
with open("result/final_state.pdbx", mode="w") as file: # pdbxはPyMOLで開けないのでcifで保存したほうがよさそう
    PDBxFile.writeFile(simulation.topology, state.getPositions(), file)

tg = time.perf_counter()
dt = tg - ts
h = dt // 3600
m = (dt % 3600) // 60
s = dt % 60
print(f"elapsed time: {h} h {m} min {s} sec")