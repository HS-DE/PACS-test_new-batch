这是我们的一个研发项目叫做PACS，不过目前我们发现这个项目有问题，结果并不好，或者说这个过程做与不做没有区别，因此我们想重新研究一下PACS流程，以下是仓库的数据信息，如果还有疑问点，可以先问我，我来做解答；

docs：存放我们的技术文档（technote）与Proteonano™ HT Quality Statement (PQS)，用于了解这个项目的情况，technote可能提到了 Inter-Study部分和Contamination Filtering部分，这在本次研究中都不做考虑；

data：包含所有用到的数据，包括蛋白定量表、样本信息、HELP的定量表；

metadata：2个表，1个是非队列样本信息all_QC_metadata.xlsx，即QC、Internal Calibrator(QC)、Blank，里面的Neat、Blank样本都可以删掉，分析用不到；另一个是队列样本信息all_sample_metadata.xlsx。样本信息中都有 sample_id、sample_name、group1/2（优先用group2吧），Plate是对应的样本所处的96孔板，Run是对应样本加样的设备批次。

pg_matrix：里面有1个蛋白定量表pg_matrix.tsv，包含非队列样本与队列样本，共244个样本，前面5列是注释列，蛋白ID主要看第一列；

在我们前期的研究中，以下样本是outlier，可以在分析中剔除这些样本：
> remove_sample
 [1] "WJ-1-5"   "WJ-1-7"   "WJ-1-19"  "WJ-1-21"  "WJ-1-44"  "WJ-1-67"  "WJ-1-90"  "WJ-1-101" "WJ-1-102"
[10] "WJ-1-103" "WJ-1-104" "WJ-1-105" "WJ-1-112" "WJ-1-146" "WJ-1-148" "WJ-1-157" "WJ-1-171" "WJ-1-176"
[19] "WJ-1-177" "WJ-1-178" "WJ-1-187" "WJ-1-203" "WJ-1-208" "WJ-1-209" "WJ-1-212" "WJ-1-213" "WJ-1-43" 
[28] "WJ-1-59"  "WJ-1-83" 

HELP：里面就是所有样本的 HELP 定量表，行名为内标肽（即HELP），列名为sample_id，数值为强度。

HELP定义：指内标肽，我们加了50条内标肽（HELP，HEavy-Labeled Peptides），最后所有样本（包括两种QC等，不包含Blank和Neat样本）共同检出的应该为25条，内标肽加的量都一样。

Internal-QC(Internal Calibrators)
- 同一个 pooled plasma (QC) 经 1 次前处理得到的肽段加入 50 PICS，在每个板上机前打 2 个重复
- 这些样本用来矫正质谱仪器状态造成的差异，即我们所说的 Instrument normalization。

QC/QC1
- 同一个 pooled plasma 样本随着每块板进行的前处理得到的肽段在上机之前加入 50 PICS
- 这个 QC 样本是用来矫正每块板之间的差异，即我们所说的 Plate normalization。

背景信息中：
Plate：板级控制（Plate Controls），每块板都会同时处理 pooled plasma QC 样本，并与研究样本并行检测分析。这些 QC 作为板级对照，用于监测孔间与板间的一致性，并用于校正由样本前处理、试剂处理以及板上操作引入的板间批次效应（plate-to-plate batch effects）。在 pooled QC 样本中同样加入 Internal Calibrators，以确保其在作为板级参考之前具备内部一致性，从而为批次对齐提供更稳健、稳定的基准。

Run：重标肽标准品（Internal Calibrators）信号检测控制（Signal Detection Controls），专用参考进样（reference injections）中加入浓度匹配的稳定同位素标记肽段作为 Internal Calibrators，用于在大量血浆/血清蛋白组研究中实现一致的定量，并具备良好的可合成可获得性。它们用于跨 Run 跟踪并校准 LC–MS 的肽段检测表现，包括信号响应、色谱稳定性以及仪器 Run-to-Run 漂移。这些 Internal Calibrators 支持信号检测层面的标准化（signal-detection normalization），用于跨 Run 的数据对齐，并通过校正与仪器相关的变异，确保整个研究期间质谱性能稳定。

script：存放脚本用；

Results：存放分析结果用；

旧示例：里面包含我们之前做过的矫正流程，给你做参考用，当时是整理好了一个RData对象，里面有所有要用到的数据，所以与现在的这种所有数据都是单独的表格文件不太一样，但是数据类型与内容形式含义都是一样的。
