#!/bin/bash

# nvidia/cuda の docker image は devel を使用する
# ダウンロードしたソースコードを解凍する
# インストール
cd gromacs-2024.2
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/gromacs_2024.2 -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DGMX_MPI=OFF -DGMX_DOUBLE=OFF -DGMX_OPENMP=ON -DGMX_GPU=CUDA
make -j16 install
source /opt/gromacs_2024.2/bin/GMXRC
gmx -version # 表示されれば成功

# リガンドトポロジーファイル作成　→　GROMACSでできないのでopenbabelとacpypeを使う
obabel -h -ipdb data/lig_4lxz.pdb -omol2 -O data/ligand.mol2

git clone https://github.com/alanwilter/acpype.git
python3 acpype/run_acpype.py -i data/ligand.mol2

# タンパク質トポロジーファイル作成
gmx pdb2gmx -f data/4lxz_receptor_protein.pdb -o data/structure.gro -ignh # amber99SB-ildnを選ぶとTIP3Pが選択可能
gmx pdb2gmx -f data/4lxz_receptor_protein.pdb -p data/structure.top 


# ファイルを統合する。詳しくは https://zenn.dev/labcode/articles/5773fe87a13e61#%E3%82%BF%E3%83%B3%E3%83%91%E3%82%AF%E8%B3%AA%E3%81%A8%E3%83%AA%E3%82%AC%E3%83%B3%E3%83%89%E8%A4%87%E5%90%88%E4%BD%93%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB%E3%81%AE%E5%90%88%E6%88%90

gmx editconf -f data/complex.gro -o data/newbox.gro -bt dodecahedron -d 1.0
gmx solvate -cp data/newbox.gro -cs spc216.gro -p data/complex.top -o data/solv.gro

gmx grompp -f data/ions.mdp -c data/solv.gro -p topol.top -o data/ions.tpr
# ions.mdp の作成。上のwebページに準拠
# トラブルシューティング　→　https://manual.gromacs.org/current/user-guide/run-time-errors.html
# forcefield の 直下に atomtypes を置く
# molecules に error で指定された差分 / 3 のSOLを追加する
gmx genion -s data/ions.tpr -o data/solvated_ions.gro -p topol.top -pname NA -nname CL -neutral

# em.mdpを作成
gmx grompp -f data/em.mdp -c data/solvated_ions.gro -p topol.top -o em.tpr
# エネルギー最小化
gmx mdrun -v -deffnm em

# リガンド-タンパク質の設定
gmx make_ndx -f data/ligand_GMX.gro -o data/index_ligand_GMX.ndx
> 0 & ! a H*
> q
gmx genrestr -f data/ligand_GMX.gro -n data/index_ligand_GMX.ndx -o data/posre_ligand_GMX.itp -fc 1000 1000 1000
> 3
# Ligand position restraints を追加
gmx make_ndx -f em.gro -o index.ndx
> 1 | 13
> q

# nvt.mdp を作成
# tc-grps を Protein_UNL SOL に変更
gmx grompp -f data/nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -deffnm nvt
# npt.mdp を作成
gmx grompp -f data/npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -deffnm npt
# md.mdp を作成
gmx grompp -f data/md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_10.tpr
gmx mdrun -deffnm md_0_10
gmx mdrun -s md_0_10.tpr -cpi md_0_10.cpt -deffnm md_0_10 # mdrunが途中終了した場合
gmx trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc mol -ur compact


# gmx_MMPBSAの構築。
# 参考は https://zenn.dev/labcode/articles/d3d109d8a7cd0e#gmx_mmpbsa%E3%81%AE%E7%92%B0%E5%A2%83%E6%A7%8B%E7%AF%89
# https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/Protein_ligand/ST/
gmx_MMPBSA --create_input gb pb ala decomp
