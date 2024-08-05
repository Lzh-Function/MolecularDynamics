# MolecularDynamics

## 参考文献

[MD計算のためのオープンソースソフトウェアをさらっと](https://magattaca.hatenablog.com/entry/2022/02/06/193011) 
[OpenMMをステップバイステップで](https://magattaca.hatenablog.com/entry/2022/02/12/201641) 
- Part10まで存在。現実行コードの基礎。
[OpenMM](https://openmm.org/) 
[openmmforcefields](https://github.com/openmm/openmmforcefields)

## インストール
Dockerfileを参照
- そのままだとnglviewがうまく回らない。可視化をすべてPyMOLに委ねるか、nglviewのセットアップに成功したcontainerを解凍して使う

### 各ディレクトリの概要
- 240612_tutorial: 参考文献を引用したタンパク質のみのMD
- 240613_tutorial_ligand: 参考文献を引用したタンパク質-リガンドドッキングのMD
- 240625_hsp90: 自ら用意したタンパク質-リガンドでのMD
- 240704_hdac2: 金属を含んだタンパク質に対するドッキング