import numpy as np
from dragon_positon import *

# 输出初始值
print(b * 2 * np.pi)  # 输出: 0.55

# 修改 b 的值
change_b(0.45 * 2 * np.pi)

# 输出修改后的值
print(b * 2 * np.pi)  # 输出: 0.45
