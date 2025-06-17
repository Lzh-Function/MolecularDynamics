# MolecularDynamics

MD用の各種スクリプト置き場
- OpenMM
- GROMACS
- Amber (pmemd, AmberTools)

### 参考文献
- **OpenMM**
    - [MD計算のためのオープンソースソフトウェアをさらっと](https://magattaca.hatenablog.com/entry/2022/02/06/193011) 
    - [OpenMMをステップバイステップで](https://magattaca.hatenablog.com/entry/2022/02/12/201641) 
        - Part10まで存在。現行コードの根幹となっている。
    - [OpenMM](https://openmm.org/) 
    - [openmmforcefields](https://github.com/openmm/openmmforcefields)
- **GROMACS**
    - [In silico創薬を論文に沿ってフォローする](https://zenn.dev/labcode/articles/b8bc3b4320dcff)
        - 現行コードの根幹となっている。
    - [GROMACS Documentation](https://manual.gromacs.org/current/index.html)
- **Amber**
    - アカデミックライセンス無償化が2023年なので易しいドキュメントは少ない...。
    - [Amber Manual](https://ambermd.org/Manuals.php)
    - [Amber Tutorial](https://ambermd.org/tutorials/)

### experimentディレクトリの概要
- `nemoto` ( [Shumpei Nemoto](https://github.com/Nemoto-S) )
    - 240612_tutorial: 参考文献を引用したタンパク質のみのMD
    - 240613_tutorial_ligand: 参考文献を引用したタンパク質-リガンドドッキングのMD
    - 240625_hsp90: 自ら用意したタンパク質-リガンドでのMD
    - 240704_hdac2: 金属を含んだタンパク質に対するドッキング
- `lzh` ( [Takuho Ri (Zehao Li)](https://github.com/Lzh-Function) )
    - 250607_amber_testrun: Amberを用いたHDAC2+vorinostat(7zzs)のMD/MMPBSA
        - `xleap`による前処理（力場適用）後にHDAC2のbinding siteが消失していたため、要再検討