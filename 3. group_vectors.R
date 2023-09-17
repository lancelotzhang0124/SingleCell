#################构建分组向量#####################
#因为我的元数据里是按Control和Suicide分好了组，所以我们按字段分组
library("stringr")
# Detect "control" or "suicide" in colnames
is_control <- str_detect(colnames(mat), "Control")
# Create a group_list based on the detection
group_list <- ifelse(is_control, "Control", "Suicide")
# Convert group_list to a factor with specified levels
group_list <- factor(group_list, levels = c("Control", "Suicide"))
group_list
