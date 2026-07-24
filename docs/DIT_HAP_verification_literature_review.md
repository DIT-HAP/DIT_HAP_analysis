# DIT-HAP 敲除验证结果的生物学意义调研

## 摘要

- **验证结果集中在少数几类保守、非冗余的"核心机器"上**：26S 蛋白酶体、血红素合成、Fe-S/CIA 组装、核糖体生物发生与翻译起始、N-/O-糖基化与细胞壁合成、囊泡运输（SNARE/TRAPPII）、剪接体 U1/tri-snRNP。这些通路在 *S. cerevisiae* 中几乎都是"教科书级必需"，而 *S. pombe* 基因组敲除文库（Kim et al. 2010, PMID:20473289；Hayles et al. 2013, PMID:23697806）却大量将其同源基因标注为 Viable——本次 DIT-HAP + 四分子解剖验证把这些错误的"WT-like/viable"标注纠正为 Essential，**是对已有认知的"纠偏"而非"颠覆"**：它让裂殖酵母的必需基因图谱与跨物种保守性重新自洽。
- **少数几个基因是真正的物种特异性发现**，而不只是修正文库错误：pre9（20S α3，budding yeast 中唯一可被 α4 替代的非必需亚基，但在裂殖酵母中必需，提示 α4-α4 补偿机制不存在或不足）、gos1/trs65（Golgi SNARE/TRAPPII 中在芽殖酵母可冗余但在裂殖酵母必需）、omh3/anp1（Golgi 甘露糖转移酶家族冗余度在裂殖酵母中似乎更低）、ifa38/css1（脂质延伸/鞘脂代谢中的物种差异）、abo1（ATAD2 家族组蛋白伴侣，其芽殖酵母同源物 Yta7 非必需）。这些是本研究**新增的生物学信息**，值得优先跟进。
- **DIT-HAP 筛选本身在一小类基因上存在系统性"假阴性"**：nse3、slx8、rmi1、usp103、usp108、rsc7、med9 这七个基因在文库/文献中早已确认必需，但 DIT-HAP 的 DR 值接近于零（0.006–0.35）。这提示 DR 读数在这些位点上未能捕捉到真实的必需性，可能与转座子插入密度、基因长度或耗竭动力学有关，而不是生物学上的例外。
- **一小部分验证结果本身与两个物种的已发表文献都冲突，应谨慎对待**：pka1（裂殖酵母 PKA 在标准条件下被广泛报道为非必需）、lsm1（细胞质 mRNA 降解因子，两个物种都非必需，且被错误归入 snRNP 剪接组）、elp1 与 fkh2（Elongator 和 forkhead 冗余在两个物种都是非必需的教科书结论）、srw1/ste9（APC/C 辅因子，经典的"可育性缺陷但存活"基因）。这些"Essential, small colonies"类调用更可能反映**严重的定量适应性缺陷**而非真正的致死性，值得用更严格的标准复核。
- 输入数据中"线粒体"组的 per-gene/claims 字段为空、"细胞周期与核心代谢单基因"组的资料被截断，这两处存在明显的**数据缺口**，下文已明确标注，不做臆测性补全。
- 就"Q0"框架（某些原认为 WT-like 的基因被验证为 Essential；某些原 Essential 的基因被验证为 WT-like/Viable）而言：本次提供的数据中**只观察到"WT-like → Essential"这一个方向**的重新分类，没有出现"Essential → WT-like/Viable"的反向案例，特此说明，不做虚构。

---

## 一、验证结果落在哪些通路/复合物，以及对复合物认识的改变（Q1）

### 1.1 概览表

| 功能模块 | 代表基因 | 重新分类方向 | 对复合物认识的改变 | 关键结论 |
|---|---|---|---|---|
| 26S 蛋白酶体（20S CP + 19S RP） | pre9, pre6, rpt4, rpt6 | WT-like→Essential | **纠偏**：恢复"整个蛋白酶体在裂殖酵母中全体必需"的预期一致性 | pre9 是唯一的**物种特异性差异**：其芽殖酵母同源物是 20S 中唯一非必需亚基（可被 α4 替代），裂殖酵母中却必需 |
| 血红素/四吡咯合成 | hem2, hem3, hem12 | WT-like→Essential | **纠偏**：消除与 hem1（已知必需）矛盾的内部不一致 | 裂殖酵母作为 petite-negative 酵母，对血红素依赖性强于可营养缺陷补救的芽殖酵母 |
| 胞质+线粒体 Fe-S/CIA 组装 | dre2, tah18, nfs1 | WT-like→Essential | **纠偏**：把 ISC（线粒体硫供体）与 CIA（胞质电子传递链）两臂的起始步骤都确立为不可绕过的单点 | 与芽殖酵母已知必需性完全对齐，无需引入新组分 |
| N-/O-糖基化、GPI、细胞壁合成 | alg14, uap1, SPAC13C5.05c, ogm2, gls1, cwh43, bgs4, omh3, anp1 | WT-like/条件依赖→Essential | **混合**：UDP-GlcNAc 核心供给纠偏；omh3/anp1/cwh43 揭示**真实物种差异**（冗余度降低） | 是本报告中"纠偏"与"新发现"并存最典型的模块 |
| SNARE/Golgi/内吞、TRAPPII | sec20, sft1, sec63, trs130, pan1, gos1, trs65 | WT-like→Essential | **混合**：5 个基因纠偏对齐芽殖酵母；gos1、trs65 是**真实物种差异** | 裂殖酵母 Golgi SNARE 束和 TRAPPII 复合物容忍冗余度更低 |
| 核糖体生物发生与翻译 | nog2, kri1, tif6, utp20, fib1, nat10, rpl14, tif213, pabp, drs1 | WT-like→Essential（nat10 已知必需） | **纠偏**：整块模块与芽殖酵母对齐 | rpl14 是"WGD 旁系同源基因掩盖必需性"的教学案例 |
| 剪接体 U1/tri-snRNP/step-1 | snu23, cwf16, usp103, usp108, (lsm1 误归类) | 混合方向 | snu23/cwf16 **纠偏**；usp103/usp108 是**DIT-HAP假阴性**（本已知必需） | lsm1 的必需性调用与两物种文献均冲突，应谨慎 |
| 转录与染色质（TFIIH/TFIIIC/RSC/Mediator/Elongator） | tfb2, sfc4, rba50, ssr3, rsc7, med9, elp1, fkh2, abo1 | 混合方向 | 核心通用转录因子**纠偏**；rsc7/med9 是**DIT-HAP假阴性**；elp1/fkh2/abo1 是**存疑的物种特异性发现** | elp1、fkh2 与两物种已知冗余性直接冲突，最需复核 |
| 脂肪酸合成/延伸、磷脂酰肌醇激酶 | fas1, elo2, ifa38, stt4, css1 | WT-like/小菌落→Essential | **纠偏+新发现**：脂质生物合成整条链在裂殖酵母中呈现"慢耗竭必需基因"特征，此前被文库/传代稀释掩盖 | ifa38、css1 是真实物种差异（芽殖酵母中非必需） |
| 信号转导（Rho1/PKA/PLC） | rho1, pka1, plc1 | WT-like/小菌落→Essential | rho1 **纠偏**（确认性）；plc1 与条件依赖性一致；**pka1 与已发表共识直接冲突** | pka1 是本报告中最值得警示的存疑案例 |
| 线粒体生物发生/mtDNA维持 | bot1, mrpl50, mrps26, mhr1, pog1, sam50, ups1, aim22, nad1, dml1 | 混合（数据不全） | Petite-negative vs petite-positive 框架解释了跨物种必需性差异 | **该组 per-gene/claims 数据缺失**，仅有组水平文献综述 |
| 基因组稳定性（HR/Smc5-6/STUbL/RecQ-Top3） | rad51, nse6, nse3, slx8, rmi1 | 混合方向 | nse3/slx8/rmi1 是**DIT-HAP假阴性**；rad51/nse6 更像**严重适应性缺陷**而非经典致死 | 展示 DR 幅度与二元必需性判定是互补而非冗余的读数 |
| 细胞周期/核心代谢单基因 | srw1, eno101, acs1, ser3, hal3, nus1, tyw3, pmt1 | 混合（数据被截断） | nus1/hal3 与已知必需复合物对齐；其余多为**存疑**（旁系同源冗余基因单独"必需"的调用需要复核） | **输入数据在该组被截断**，claims 缺失 |

### 1.2 最强的几个故事详解

**26S 蛋白酶体（pre9, pre6, rpt4, rpt6）**：20S 核心颗粒的两个外圈 α 亚基（pre9=α3, pre6=α4）与 19S 调节颗粒的两个 AAA+ ATPase（rpt4, rpt6）此前在 PomBase 上分别被标注为 Viable / Viable / Viable / "Depends on conditions"（source: https://www.pombase.org/gene/SPAC13C5.01c ; https://www.pombase.org/gene/SPBC106.16 ; https://www.pombase.org/gene/SPCC1682.16 ; https://www.pombase.org/gene/SPBC23G7.12c）。但蛋白酶体在真核生物中高度保守且几乎全体必需——芽殖酵母中除 α3/PRE9 外所有亚基必需，α3 之所以是例外，是因为 α4/Pre6 可以替代它组成 α4-α4 蛋白酶体（Velichutina et al. 2004, source: https://pubmed.ncbi.nlm.nih.gov/14739934/）。rpt6 本身在裂殖酵母中早已是经典致死基因 let1（Michael 1994, source: https://pubmed.ncbi.nlm.nih.gov/8056332/），pre6 携带 mts7 同义名并被 Penney et al. 2012 证实为必需（source: https://core.ac.uk/download/pdf/28975224.pdf）。DIT-HAP+四分子将全部四个基因确认为 Essential，**修复了文库的假阴性标注，让蛋白酶体重新读作一个整体必需的复合物**。最有信息量的是 pre9：它的芽殖酵母同源物是唯一非必需的 20S 亚基（有 α4 补偿），而在裂殖酵母中却必需——提示这种补偿机制在裂殖酵母中不存在或不足以维持存活，这是一个具体、可检验的物种差异，而不仅仅是文库伪影。

**血红素合成（hem2, hem3, hem12）**：这是一条线性、无分支的八步通路（hem1→hem2→hem3→ups1→hem12→hem13→hem14→hem15），此前 hem1 已被证实为必需基因（Normant et al. 2018, source: https://pmc.ncbi.nlm.nih.gov/articles/PMC5925805/），而 hem2/hem3/hem12 却被文库标为 Viable——这在逻辑上自相矛盾（阻断线性通路任一步骤都应该同样致死）。芽殖酵母中 HEM2、HEM3 因血红素营养缺陷而必需，HEM12 则因为培养基可补救而"非必需但生长极慢"（source: https://www.yeastgenome.org/locus/S000003008 ; https://www.yeastgenome.org/locus/S000002364/phenotype ; https://www.yeastgenome.org/locus/S000002454）。裂殖酵母是 petite-negative 酵母，缺乏芽殖酵母那种可补救的营养缺陷路径（source: https://pmc.ncbi.nlm.nih.gov/articles/PMC10582590/），因此三者被验证为 Essential 恰好**消除了此前的内部矛盾**，确立整条血红素通路在裂殖酵母中全程必需。

**Fe-S/CIA 组装（dre2, tah18, nfs1）**：这是本组中"文库标注与跨物种预期最直接冲突"的案例——三者在芽殖酵母中都是必需基因（Nfs1 是硫供体去硫酶，Tah18-Dre2 是启动胞质 Fe-S 组装的电子传递链，Vernis et al. 2009 称其为"新发现的必需复合物"，source: https://www.pombase.org/reference/PMID:19194512），而 PomBase 却把三个裂殖酵母同源物全部标为 Viable（source: https://www.pombase.org/gene/SPBC337.10c ; https://www.pombase.org/gene/SPAC1296.06 ; https://www.pombase.org/gene/SPBC21D10.11c）。验证结果把三者都纠正为 Essential，说明 ISC（线粒体）与 CIA（胞质）两条臂的起始步骤在裂殖酵母中同样是不可绕过的单点故障，不存在跨物种的例外冗余。

**N-/O-糖基化、GPI、细胞壁（9 个基因）**：这组同时包含"纠偏"和"新发现"两类故事。alg14、uap1、SPAC13C5.05c 三者是 UDP-GlcNAc 供给与 N-糖基化起始的核心酶，其芽殖酵母同源物（ALG13/ALG14、QRI1、PCM1/AGM1）均明确必需（source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2583287/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC2747872/），文库的 Viable 标注是明显的假阴性。bgs4、ogm2、gls1 此前文献/文库已知条件依赖或必需（bgs4 在胞质分裂中必需，Cortés et al. 2005, source: https://journals.biologists.com/jcs/article/118/1/157/28237/），验证只是确认。而 **omh3（Golgi α-1,2-甘露糖转移酶，KRE2/KTR 家族）、anp1（Golgi mannan polymerase I 亚基）、cwh43（GPI 脂质重塑酶）** 在芽殖酵母中都因旁系同源基因冗余或本身非必需而可存活，却在裂殖酵母被验证为必需——这提示**裂殖酵母 Golgi 糖基转移酶家族的功能冗余度低于芽殖酵母**，是一个具体、值得后续做旁系同源基因互补实验验证的假说。

**SNARE/Golgi/TRAPPII/内吞（7 个基因）**：sec20、sft1、sec63、trs130、pan1 五者的芽殖酵母同源物（SEC20、SFT1、SEC63、TRS130、PAN1）全部必需（source: https://www.yeastgenome.org/locus/YDR498C ; https://www.yeastgenome.org/locus/YOR254C ; https://www.yeastgenome.org/locus/YMR218C ; https://www.yeastgenome.org/locus/YIR006C），但裂殖酵母文库全部标为 Viable——验证纠正了这五处假阴性。真正有信息量的是 **gos1 和 trs65**：Gos1 在芽殖酵母中明确非必需（McNew et al. 1998，"GOS1 was not an essential gene"，source: https://pubmed.ncbi.nlm.nih.gov/9755865/），Trs65 也非必需（source: https://www.yeastgenome.org/locus/YGR166W）；若验证的必需性成立，说明裂殖酵母的 Golgi SNARE 束（Sed5-Ykt6-Gos1-Sft1）与 TRAPPII 复合物比芽殖酵母**容忍更少的亚基冗余**，是真实的物种分化而非标注错误（trs65 本身在文库中已经是 Inviable，所以这里验证只是确认）。

**核糖体生物发生与翻译（10 个基因）**：这是最"整齐"的纠偏案例——nog2、kri1、tif6、utp20、fib1、nat10、tif213、pabp、drs1 的芽殖酵母同源物全部必需（NOG2、KRI1、TIF6、UTP20、NOP1、KRE33、GCD11、PAB1、DPS1，source 分别见 https://www.yeastgenome.org/locus/YNR053C 等），文库却给出 9 个 Viable 标注，只有 nat10 本已是 Inviable。**rpl14 是本组最具教学价值的案例**：芽殖酵母中 L14 由两个全基因组复制（WGD）旁系同源基因 RPL14A/RPL14B 编码，单基因缺失可存活、双缺失致死（source: https://www.yeastgenome.org/locus/YHL001W），这种冗余掩盖了该基因家族本质上的必需性；裂殖酵母的 rpl14 若为单拷贝功能基因，其必需性调用恰好揭示了"WGD 冗余如何系统性地压低必需基因清单的表观覆盖率"。

**剪接体 U1/tri-snRNP（snu23, cwf16, usp103, usp108；lsm1 误归类）**：snu23（tri-snRNP）与 cwf16（NTC 相关 first-step 因子）的芽殖酵母同源物 SNU23、YJU2/CWC16 均必需（source: https://www.yeastgenome.org/locus/YDL098C ; https://www.yeastgenome.org/locus/YKL095W），而裂殖酵母文库标为 Viable（cwf16Δ 甚至被 Sasaki-Haraguchi et al. 2015 作为可存活株研究，source: https://www.pombase.org/reference/PMID:26302002）——验证把两者纠正为 Essential。usp103（U1-C，保守亚基）与 usp108（裂殖酵母特有 U1 蛋白）**已经**在 PomBase 中被标为 Inviable（Newo et al. 2007, source: https://www.pombase.org/reference/PMID:17264129），但它们在 DIT-HAP 筛选中的 DR 极低（0.155、0.273）——说明这两处不是文库错误，而是**筛选本身的假阴性**。lsm1 是一个明显的分组错误：它是细胞质 mRNA 降解（P-body/Lsm1-7-Pat1）因子而非剪接体亚基，在两个物种中都非必需（source: https://www.yeastgenome.org/locus/YJL124C），任何把它判为必需的结论都需要独立复核。

**转录与染色质（9 个基因）**：tfb2（TFIIH）、sfc4（TFIIIC）、rba50（Pol II 组装因子）、ssr3（RSC/SWI-SNF 共享亚基）的芽殖酵母同源物均必需，文库标为 Viable——验证纠正。rsc7、med9 文库中**已经**是 Inviable，但 DR 极低（0.345、0.092），是**DIT-HAP假阴性**；有意思的是它们各自的芽殖酵母同源物（Npl6、Cse2）却都非必需（source: https://www.yeastgenome.org/locus/S000005293），说明这两个基因是"裂殖酵母中真正必需但芽殖酵母中非必需"的物种差异，DIT-HAP 只是没能测出其真实的必需性幅度。**elp1（Elongator）、fkh2（forkhead）和 abo1（ATAD2 家族组蛋白伴侣）是最需要谨慎对待的三个基因**：Elongator 在两个酵母中都是公认非必需的（source 见 elp1 词条列表），forkhead 家族在芽殖酵母中因 Fkh1/Fkh2 冗余而非必需，Yta7（abo1 同源物）在芽殖酵母中也非必需——三者若真的在裂殖酵母中必需，将是意料之外的新发现，但同样可能是验证方法本身（菌株背景、培养条件、小菌落误判）导致的假阳性，原始数据本身也明确将 elp1/fkh2 标注为"最需要复核"。

**脂肪酸合成/延伸与磷脂酰肌醇激酶（fas1, elo2, ifa38, stt4, css1）**：这组的叙事是"慢耗竭必需基因被稀释法筛选/文库掩盖"——脂质合成缺陷细胞往往依靠继承的脂质库和培养基中可清除的脂肪酸维持数代分裂才表现出致死性，这正是文库单次传代观察容易漏检、而 DIT-HAP 的多代竭尽读数（加上四分子解剖）能够揭示的场景。stt4 在两个物种中都必需，是最强确认性证据（source: https://www.yeastgenome.org/locus/YLR305C ; PomBase 已标注 stt4Δ inviable）。ifa38 和 css1 在芽殖酵母中非必需（source: https://www.yeastgenome.org/locus/YBR159W ; https://www.yeastgenome.org/locus/YER019W），若验证的必需性成立，则是真实的物种差异，提示裂殖酵母对脂肪酸延伸循环和鞘脂降解的依赖性高于芽殖酵母。

**信号转导（rho1, pka1, plc1）**：rho1 作为 (1,3)-β-葡聚糖合酶复合物的调节亚基，本身在两个物种中都已知必需（source: https://www.yeastgenome.org/locus/S000006369 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC6802929/），验证只是确认性的。plc1 与已知的条件依赖表型（viable-but-sick 到 inviable population）一致。**pka1 是这组、也是本报告中最值得标记的矛盾点**：多篇文献明确指出裂殖酵母 PKA 在标准培养条件下非必需（source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2941774/），PomBase 也标注为 Viable（source: https://www.pombase.org/gene/SPBC106.10），而验证却给出 Essential——这与已发表共识直接冲突，应被视为存疑结果，而非确定的新发现（详见第三节）。

**线粒体生物发生/mtDNA维持（10 个基因）**：**该组的 per_gene 和 claims 字段在输入数据中为空**，只有组水平的 essentiality_consensus 文本可用。该文本给出的核心框架是：裂殖酵母是 petite-negative 酵母（不能在缺失 mtDNA 的情况下存活），这与 petite-positive 的芽殖酵母形成鲜明对比，解释了为什么 sam50（Kozjak et al. 2003, PMID:14570913）、dml1（PMID:12702300）、nad1/FAD1（PMC231949）在芽殖酵母中必需，而线粒体核糖体亚基、pol-γ（pog1/MIP1）、脂酰化连接酶（aim22/AIM22-LIP3）在芽殖酵母中却非必需（因为可以靠发酵存活）。同时该文本指出裂殖酵母中许多线粒体基因是"低糖条件下条件必需"（PMID:37859837），这是一个重要的补充维度，但**由于缺少逐基因的 DR 值和验证结果，无法在此对单个基因下确定性结论**，建议在后续分析中补齐这部分数据。

**基因组稳定性（rad51, nse6, nse3, slx8, rmi1）**：nse3（Smc5-6 核心亚基）、slx8（STUbL 催化亚基）、rmi1（RecQ-Top3-Rmi1 复合物）在裂殖酵母文库中**已经**是 Inviable，但 DR 几乎为零（0.006、0.156、0.006）——这是清晰的**DIT-HAP假阴性**，说明筛选在这三个位点系统性漏检了已知的必需性。相比之下，rad51（HR 核心重组酶）和 nse6（Smc5-6 辅助亚基 Nse5-Nse6 二聚体）在两个物种的经典遗传学文献中都是**非必需**（source: https://www.yeastgenome.org/locus/YER095W ; Pebernard et al. 2006 明确指出 Nse5/Nse6 在裂殖酵母中非必需，source: https://pmc.ncbi.nlm.nih.gov/articles/PMC1430260/），验证给出的"Essential, small colonies"更可能反映的是**严重的定量适应性缺陷**而非经典的二元致死性——这提醒我们四分子解剖+菌落面积法在评分"小菌落"时，本质上捕捉的是一个连续的适应性谱，把它折叠成二元的 Essential 标签可能夸大了结论的强度。

**细胞周期与核心代谢单基因（srw1, eno101, acs1, ser3, hal3, nus1, tyw3, pmt1）**：**该组的输入数据在文本中被截断**（在描述 srw1/ste9 的矛盾之处时中断），claims 数组未提供。可以确认的是：nus1（DDS 顺式异戊烯基转移酶复合物核心亚基）和 hal3（CoA 生物合成 PPC 脱羧酶异源三聚体亚基）与芽殖酵母中已知必需的复合物（source 中提及 UniProt Q12063；Olzhausen et al. 2013, PMID:23789928）逻辑一致，属于"纠偏"。而 eno101（enolase, 有旁系同源基因 eno102）、ser3（有旁系同源基因，双缺失才是丝氨酸营养缺陷）、acs1（其旁系同源基因 acs2 在葡萄糖存在时必需，双缺失才致死，source: PMID:7649171）、srw1/ste9（APC/C 辅活化因子，经典可育但存活基因）在芽殖酵母中都是"因旁系同源冗余而单独非必需"的基因——若这些基因在裂殖酵母中被判定为单独必需，需要格外谨慎，因为这意味着该基因是功能上占主导的拷贝，这是比"纠正文库标注错误"更强的断言，需要更直接的实验证据（如互补实验）支持。**由于该组的具体证据链未在输入数据中完整给出，这部分结论应视为初步方向，而非已验证的结论。**

---

## 二、功能未被充分研究的基因及跨物种同源物信息（Q2）

以下基因在 *S. pombe* 中此前研究较少（无正式基因名、无同源物、或被归入更大家族而缺乏个体研究），本次验证结合跨物种同源物功能信息，能够显著提升对它们的认识：

**SPAC13C5.05c（N-acetylglucosamine-phosphate mutase）**：这是一个至今没有正式基因名的基因，本身在裂殖酵母文献中几乎没有独立研究。但其芽殖酵母同源物 PCM1/AGM1（YEL058W）被 PomBase 直接标注为"Essential"（source: https://www.pombase.org/gene/SPAC13C5.05c），人类同源物是 PGM3（先天性糖基化障碍相关基因）。这是 UDP-GlcNAc 生物合成通路（N-糖基化、GPI 锚定、几丁质合成的共同前体供给）中的关键酶。此前文库把它标为 Viable，是明显的假阴性。**我们的验证（Essential）+ 跨物种同源物信息，把一个"无名"基因直接升级为一个功能明确的、保守的必需代谢酶**，值得考虑给它命名并纳入 UDP-GlcNAc 通路的核心必需基因清单。

**usp108（U1 snRNP 相关蛋白，裂殖酵母特有）**：这是裂殖酵母特有、无芽殖酵母同源物的 U1 snRNP 组分，最初由 Newo et al. 2007 的蛋白质组学研究鉴定为"三个必需的、物种特有的 U1 蛋白"之一（source: https://www.pombase.org/reference/PMID:17264129）。它在 PomBase 中已经是 Inviable，但 DIT-HAP DR 只有 0.273——本次验证的价值不在于重新分类，而在于**用独立的遗传学方法（四分子解剖）确认了 2007 年蛋白质组学推断的必需性**，同时也暴露出 DIT-HAP 筛选在这个位点的检测盲区，为理解裂殖酵母剪接体的"物种特化附件"提供了双重证据支撑。usp103（保守的 U1-C/YHC1 同源物）情况类似，是保守核心亚基的确认而非新发现。

**omh3、anp1（Golgi α-甘露糖转移酶家族）**：这两个基因分别是 KRE2/KTR 家族和 mannan polymerase I 复合物的成员，在芽殖酵母中因为家族内部高度冗余而个体非必需（source: https://www.pombase.org/gene/SPCC777.07 ; https://www.pombase.org/gene/SPBC1734.04）。它们在裂殖酵母中此前也几乎没有单独的功能研究。如果验证的必需性成立，**这是一个具体、可检验的假说：裂殖酵母的 Golgi 甘露糖转移酶家族成员数量更少或功能重叠度更低**，值得通过检索裂殖酵母基因组中该家族的旁系同源基因数目来验证（若裂殖酵母该家族确实比芽殖酵母缩小，就为必需性提供了直接的机制解释）。

**cwh43（GPI 脂质重塑酶）**：芽殖酵母和人类的同源物都非必需（source: https://www.pombase.org/gene/SPAC589.12），此前 PomBase 记录了 cwh43 缺失后的"swollen-elongated cell / abolished-cytokinesis"表型但整体归为 Viable。若验证的必需性成立，说明**裂殖酵母对 GPI 锚的神经酰胺重塑步骤的依赖程度超过芽殖酵母**，可能与裂殖酵母细胞壁结构（β-葡聚糖为主，甘露聚糖层结构不同）对 GPI 锚定蛋白的稳定性要求更高有关。

**abo1（ATAD2 家族组蛋白伴侣）**：其芽殖酵母同源物 Yta7（YGR270W）明确非必需（source: https://www.yeastgenome.org/locus/S000003502），人类同源物 ATAD2/ATAD2B 是肿瘤生物学中被广泛研究的溴结构域蛋白（癌症药物靶点）。此前裂殖酵母中对 abo1 的研究非常有限。**如果本次的必需性调用能被独立复核确认，这将是本次分析中最具新颖性的发现之一**——它意味着裂殖酵母把一个在人类肿瘤研究中重要、但在酵母遗传学中被视为"非必需伴侣"的因子，变成了不可或缺的组蛋白伴侣，这个物种特异性差异本身值得深入研究其分子机制（例如裂殖酵母是否缺少其他冗余的 H3-H4 伴侣）。

**gos1、trs65（Golgi SNARE / TRAPPII 附属亚基）**：两者在芽殖酵母中都明确非必需（McNew et al. 1998 对 gos1 的经典研究，source: https://pubmed.ncbi.nlm.nih.gov/9755865/；trs65 见 https://www.yeastgenome.org/locus/YGR166W），此前在裂殖酵母中也缺乏专门研究（trs65 除外，它在文库中本已是 Inviable）。**验证结果为 gos1 提供了新的必需性信息**，提示裂殖酵母 Golgi SNARE 复合物（Sed5-Ykt6-Gos1-Sft1）中 Gos1 承担的角色可能比其芽殖酵母对应物更加不可替代，是值得跟进的物种分化案例。

**rpl14**：芽殖酵母中因 WGD 旁系同源基因 RPL14A/RPL14B 的冗余而在遗传学分析中"隐身"（单基因删除不致死），这使得该基因家族的必需性长期被低估。裂殖酵母中关于 rpl14 单独功能的文献报道很少。本次验证为其必需性提供了直接证据，也是理解"基因组复制事件如何系统性地压低表观必需基因清单覆盖率"的一个具体例证。

**pre9（proteasome α3）**：虽然其芽殖酵母同源物研究充分，pre9 在裂殖酵母中的必需性此前缺乏直接遗传学验证。**关键的新认识是：芽殖酵母中允许 α4-α4 蛋白酶体绕过 α3 缺失（Velichutina et al. 2004），而裂殖酵母中这一补偿机制显然不起作用或不充分**——这为后续研究"裂殖酵母蛋白酶体组装是否存在这种柔性替代机制"提供了明确的实验切入点（例如检测 pre9 缺失细胞中是否存在异常的蛋白酶体组装中间体）。

**tyw3、pmt1（tRNA 修饰酶）**：pmt1 是 Dnmt2/TRDMT1 家族的 tRNA C5-胞嘧啶甲基转移酶，在芽殖酵母中**没有同源物**，其人类同源物 TRDMT1 在应激反应和肿瘤生物学中有一定研究基础；tyw3（wybutosine 生物合成）在芽殖酵母中非必需（生长缓慢）。这两个基因在裂殖酵母中此前研究都很有限，**由于本组输入数据被截断，缺乏具体的 claims/URL 支撑**，此处仅作为方向性提示，不作确定性结论。

---

## 三、验证失败/不一致基因的原因归因（Q3）

将本数据集中所有"意外"或"不一致"的验证结果，按数据 vs 生物学原因分为三类：

### (A) 文库假阴性（deletion-library false negatives）——占主导地位

这是数据集中**数量最多、证据最扎实**的一类：文库（Kim et al. 2010 / Hayles et al. 2013）把大量核心保守基因错误标注为 Viable，而其芽殖酵母同源物早已明确必需。属于这一类的基因包括：pre9、pre6、rpt4、rpt6（蛋白酶体）；hem2、hem3、hem12（血红素合成）；dre2、tah18、nfs1（Fe-S/CIA）；alg14、uap1、SPAC13C5.05c（糖基化起始）；sec20、sft1、sec63、trs130、pan1（囊泡运输）；nog2、kri1、tif6、utp20、fib1、tif213、pabp、drs1、rpl14（核糖体生物发生/翻译）；snu23、cwf16（剪接）；tfb2、sfc4、rba50、ssr3（转录）；fas1（脂肪酸合成）。**这些不是"验证失败"，而是验证成功纠正了文库的系统性错误**。可能的技术原因包括：竞争性池式筛选中的营养交叉喂养（例如血红素/ALA 交叉营养）、慢生长突变体被误判为存活、条码错配或不完全缺失导致的假阴性。这类错误的规模之大提示，**裂殖酵母基因组敲除文库对"高度连通、保守必需复合物核心亚基"存在系统性的检测盲区**，这本身是一个值得单独报告的元结论。

### (B) DIT-HAP 筛选假阴性（screen false negatives）——一个明确但较小的类别

nse3、slx8、rmi1、usp103、usp108、rsc7、med9 这七个基因在文库/文献中**已经**被正确标注为必需（Inviable），但 DIT-HAP 的 DR 值极低（0.006–0.35），几乎与非必需基因无法区分。这明显是**DIT-HAP 筛选方法本身的技术局限**，而非文库或生物学的问题。可能的原因（数据中未直接给出，属于推测但方向合理）：转座子/敲低插入位点密度不足、基因长度过短导致插入机会少、蛋白质耗竭动力学过慢导致表型在筛选窗口内未充分显现、或者这些基因产物的功能丧失需要更长时间才能转化为生长劣势。**这类基因提示 DR 数值本身不能单独作为必需性的排除依据**，需要与独立的遗传学验证（如本次的四分子解剖）联合使用。

### (C) 与已发表文献直接冲突、需要谨慎复核的验证结果——真正意义上"验证结果本身存疑"的一类

这是最接近字面意义上"验证失败"的类别——不是文库错、也不是筛选漏检，而是**四分子解剖/菌落面积法本身给出的必需性判定，与两个物种现有的坚实遗传学文献相矛盾**：

- **pka1**：多篇独立文献明确指出裂殖酵母 PKA 在标准培养条件下非必需（source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2941774/），PomBase 标注为 Viable（source: https://www.pombase.org/gene/SPBC106.10），芽殖酵母中三个旁系同源基因需要三重缺失才致死。验证给出 Essential 与这一共识直接冲突，最可能的解释是**菌株背景差异、培养条件（如特定应激/营养限制条件下 PKA 缺失细胞的生长优势消失）或四分子解剖过程中对"极小菌落"的判定阈值设置**——不建议在没有独立重复实验的情况下把 pka1 视为确证的必需基因。
- **lsm1**：细胞质 mRNA 降解因子，在两个物种中都明确非必需，且从功能上讲根本不属于剪接体 snRNP（被误归类）。任何将其判定为必需的结果都应首先怀疑是**株系/菌落评分误差或基因型确认问题**，而非真实生物学。
- **elp1、fkh2**：Elongator 复合物和 forkhead 转录因子家族的非必需性/冗余性在芽殖酵母遗传学中是最确立的结论之一。若裂殖酵母中这两个基因确实必需，这将是显著的物种分化，但由于与已知冗余机制（fkh1/fkh2 冗余）的直接矛盾，**在原始数据中已被明确标记为"最需复核"**，应优先安排独立的敲除/回补实验验证，而非直接采纳为结论。
- **rad51、nse6**：均为"Essential, small colonies"，且两个物种的经典遗传学都认为它们非必需（缺失细胞存活但对 DNA 损伤敏感）。这类结果更准确的解读是**严重的定量适应性缺陷，而非二元的致死性**——把菌落大小的连续谱强行折叠为二元 Essential/Viable 标签，本身就是一种方法学层面的信息损失，而不是"数据错误"，需要在后续报告中用适应性评分（而不是二元标签）来呈现这类结果。

### 总体归因结论

在能够明确归因的案例中，**(A) 文库假阴性占绝大多数**，是本次验证分析最核心、最可靠的价值所在——它系统性纠正了大量保守必需复合物核心亚基在裂殖酵母文库中的错误标注。**(B) DIT-HAP 筛选假阴性**是一个规模较小但真实存在的技术局限，集中体现在几个特定基因家族（Smc5-6、STUbL、RecQ-Top3、剪接体附属因子），提示筛选灵敏度在这些位点不足，DR 数值需要谨慎解读，不能单独作为排除必需性的依据。**(C) 与已知文献直接冲突的存疑结果**数量最少，但风险最高——这些案例最可能反映的是验证方法本身（菌落评分标准、培养条件、菌株背景）的局限而非新生物学，在对外报告或发表这些"新发现"之前，应作为最高优先级安排独立复核实验。

---

## 四、总体结论与对该验证分析价值的评估

本次 DIT-HAP + 四分子解剖验证的**最大价值**在于系统性纠正了 *S. pombe* 基因组敲除文库（Kim 2010 / Hayles 2013）对大量高度保守、通常在芽殖酵母中必需的复合物核心亚基的错误"Viable"标注——涵盖蛋白酶体、血红素合成、Fe-S/CIA 组装、核糖体生物发生、囊泡运输、剪接体核心组分等几乎所有真核生物细胞的"基础设施"通路。这类纠正的证据链最扎实：几乎每一个案例都有明确的芽殖酵母同源物必需性文献支持，跨物种一致性本身就是强有力的验证。

**第二层价值**在于揭示了若干真实的物种特异性差异（pre9 的 α4 补偿缺失、gos1/trs65 的 SNARE/TRAPPII 冗余度降低、omh3/anp1/cwh43 的 Golgi 糖基化冗余度降低、abo1 的组蛋白伴侣必需性），这些是本次分析真正产生的新生物学假设，值得优先安排后续实验（如互补实验、旁系同源基因家族大小比对）来巩固。

**第三层，也是需要研究者特别注意的**，是筛选方法本身的局限性：DIT-HAP 在少数特定基因家族上存在系统性假阴性（Smc5-6/STUbL/RecQ-Top3/剪接体附属因子），而四分子解剖+菌落面积法在评分严重适应性缺陷（而非经典致死）时可能夸大结论强度，且在极少数案例（pka1、lsm1、elp1、fkh2）中给出了与已发表文献直接冲突的结果，这些应被视为需要独立重复验证的存疑发现，而非最终结论。

整体而言，这次验证分析的科学价值是明确且可防御的：它把裂殖酵母的必需基因图谱在多个保守通路上重新与跨物种生物学对齐，同时诚实地暴露了自身方法的边界所在。

---

## 附：证据强度与被驳斥/存疑的说法

以下说法应被明确标注为**低置信度、与已发表文献直接冲突、或需要独立复核**，读者不应在未经进一步验证的情况下将其视为确证结论：

1. **pka1 判定为 Essential**：与 PomBase 标注（Viable, source: https://www.pombase.org/gene/SPBC106.10）及多篇文献（source: https://pmc.ncbi.nlm.nih.gov/articles/PMC2941774/）中"裂殖酵母 PKA 在标准条件下非必需"的共识直接冲突。**不建议采信为确定的新发现**，需要独立重复实验。
2. **lsm1 若被判定为 Essential**：与两个物种（source: https://www.yeastgenome.org/locus/YJL124C）的已知非必需性均冲突，且该基因在功能上被误归入剪接体 snRNP 分组（它实际上是细胞质 mRNA 降解因子）。**应视为分组或验证过程中的潜在误判**。
3. **elp1、fkh2 判定为 Essential**：与芽殖酵母中 Elongator 复合物和 forkhead 家族（Fkh1/Fkh2 冗余）非必需性的教科书级结论冲突。原始数据本身已标记这两处"最需复核"，**在独立验证前不应作为确定结论引用**。
4. **rad51、nse6、pka1、css1 的"Essential, small colonies"判定**：这些标签把连续的适应性缺陷谱折叠为二元必需性标签，容易夸大结论强度。**建议在后续报告中改用定量适应性评分而非二元 Essential/Viable 标签**。
5. **omh3、anp1、cwh43、gos1、ifa38 的物种特异性必需性**：虽然生物学上合理（芽殖酵母中存在旁系同源冗余或本身非必需），但目前**仅有本次验证这一组数据支持**，缺乏独立的第二方证据（如互补实验、旁系同源基因缺失表型的直接比较），应标注为"初步假设，有待巩固"而非确证结论。
6. **"细胞周期与核心代谢单基因"组（srw1, eno101, acs1, ser3, hal3, nus1, tyw3, pmt1）的多数结论**：由于输入数据在描述该组时被截断，缺少完整的 claims/URL 支撑，**该组除 nus1、hal3 外的结论均应视为方向性提示，不作为已核实的结论使用**。
7. **"线粒体生物发生"组的逐基因结论**：由于输入数据的 per_gene 和 claims 字段为空，**该组目前只能在通路/物种层面（petite-negative vs petite-positive）讨论，无法对 bot1、mrpl50、mrps26、mhr1、pog1、sam50、ups1、aim22、nad1、dml1 中任何单个基因给出确定性的验证结论**，需要补充数据后重新评估。
