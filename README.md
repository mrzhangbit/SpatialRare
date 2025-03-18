# SECIA
# SECIA Pipeline Usage Guide (English & 中文)

This repository includes scripts for data analysis and consensus assessment. Please follow the steps outlined below to use the pipeline correctly.

本仓库包含用于数据分析和一致性评估的脚本，请按照以下步骤正确使用。

## Environment Requirements (环境要求)

- R version 4.3.3 or higher. 
- Python version 3.8 or higher 

## Library Dependencies (所需库)

### Python Libraries

#### For `secia-entropy.py`:
```python
import pandas as pd
import numpy as np
import os
```

#### For `secia-consensus.py`:
```python
import pandas as pd
import numpy as np
```

### R Libraries

#### For `secia.r`:
```R
library(Matrix)
library(Seurat)
library(irlba)
library(igraph)
library(rflann)
library(e1071)
library(RANN)
library(stats)
```

#### For `secia-cci.r`:
```R
library(cluster)
library(Seurat)
library(dplyr)
library(ggplot2)
library(CellChat)
library(ggalluvial)
library(svglite)
library(scran)
library(DT)
library(Matrix)
```

## Step-by-Step Instructions (运行步骤)

### Step 1: Run `secia-cci.r` 

### Step 2: Run `secia-entropy.py` 

### Step 3: Run `consensus.py` 

Place the files generated in the previous steps (`cci.csv` and `entropy.csv`) into the current directory and run:
This script performs a consensus assessment and outputs the results as consensus.csv.

### Step 4: Final Analysis  Run `SECIA.r` 


## Important Notes (注意事项)

- Ensure each step generates the required files correctly, and the files are placed in the appropriate directory.
- 确保每一步都正确生成所需文件，并将文件放在对应目录中。
- Check library dependencies and their installation status if errors occur.
- 若出现错误，请检查库的依赖关系和安装情况。

## Contact (联系)

If you encounter any issues, please open an issue on GitHub.

如遇问题，请在GitHub上提出 issue。

