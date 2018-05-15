# classify
### aims  
update theme map  

### infomation of code 
all is the interface which include these should be complied by IDL
this process can work in IDL8.3 with opening of ENVI.
### overall idea  
1. 基于待更新影像不会发生颠覆性的变化的假设（Chinese）
2. 系统地考虑了影像自身聚类类别信息，并有效的结合已有专题图对中间结果进行控制（Chinese）
1. we suppose there are overturning change between contiguous theme maps
2. the idea consider spatial distribution of classes in image sysmaticlly. And we  has a control over intermedia processes.

# TODO  
    1. we don't prepare a interface for this alogorithm but a procedure of overall flows, which is "all.pro"  
    There are main stream of this function. If you like, you can find branches easily.  
    
 
    2. we want add more explanation for this code in the future. There are some explanation in chinese now neither detailed nor sysmatic.
    we hope this doesn't make you feel trouble too much.


# which situation does it suitable
this work is prepare for wetlands mapping updating.  
It's a very dynamic object, so this work may be a good option when those need pick sample carefully withou labor involving in it.

# experience data  
some available data will be presented laterly.

# FAQ
    Q: you take into consideration of metric of auto-correction, but how does it work
    A: we suppose there are overturning change between contiguous theme maps. Therefore, we can build certain connection between same or similar classes over different time for their spatial distribution. 
    What's more, the metric is carried out under condition of whole image. In other word, the noise of correction bring by change can be smoothed in a greater scale, such as a whole image span of Landsat.
    complementary, the general metric method is similarity of spectrum. But it brings some errors for monitering wetland for its high dynamic change.
    
    
