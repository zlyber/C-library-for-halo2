# Readme

**SDU    张芷源    2023/7/9**

​        **本工程实现了Halo2中NTT、MSM以及SETUP函数的C版本并在原工程中对接。**

​        **halo2文件夹下是本工程全部代码，已经将原本的计算函数替换成了C接口版本；**

​		**MSM文件夹下是C-MSM代码**，其中main函数的输入是我截取的某次原程序（K=5）中的输入，用于示例；

​		**NTT文件夹下是C-NTT代码**，其中main函数的$\omega$是我截取的原程序（K=5）中的值，用于示例。



### 测试说明

​		**测试环境：ubuntu 18.04/20.04+vscode+rust-analyzer**

​		**测试程序：halo2/halo2_proofs/tests/plonk_api.rs** 

​		**测试方法：**如果使用vscode，先安装rust-analyzer扩展，等待扩展激活后，点击程序中的test/debug即可进行测试或debug，首次运行时，需要等待一定时							间进行编译

​		**如何更改规模：**在测试程序中的**plonk_api函数内**，可以直接更改变量K，n=2^K，默认K=5，若更改后测试不能pass，是因为没有更改对应的原根，目前我存入了K=5、10、15、20--24的原根，可以在halo2_proofs/src/arithmetic.rs中的best_fft函数内手动更改，也可以不用管

​		**C程序位置：**在halo2/halo2_proofs/lib文件夹中，hello.c存放所有C相关代码，hello.h是NTT、MSM和setup的函数声明

​		==这个版本的NTT接口已经优化，几乎没有数据类型转换时间，MSM接口未优化==

