#! /usr/bin/env python

from spdkit import gosh
import os
import time

# bbm 目录全路径
bbm_dir = os.path.abspath("vasp-opt-bbm")
# bbm_dir = os.path.abspath("sp")
cif_root_dir = os.path.abspath("jobs")

# 在 cif 文件所在的目录, 根据 cif 文件名自动生成计算目录, 比如
# 200.cif, 将在 200 目录下执行如下 bbm 命令
# bbm -t /path/to/vasp-opt-bbm 123.cif
commands = []
for f in os.listdir(cif_root_dir):
    jobname, extension = os.path.splitext(f)
    if extension in [".cif"]:
        work_dir = os.path.join(cif_root_dir, jobname)
        # os.makedirs(work_dir, exist_ok=True)
        cif = os.path.join(cif_root_dir, f)
        cmd = f"bbm -t {bbm_dir} {cif}"
        commands.append((cmd, work_dir))

# 构造任务池, 提交到远程计算调度中心, 分发到多个计算节点并行计算
#
# scheduler 服务地址默认从当前目录下 gosh-remote-scheduler.lock 文件中读取.
# 该文件由 gosh-remote 启动 scheduler 自动生成.

# work around NFS sync issue for newly created work directories
# time.sleep(5)
hub = gosh.JobHub.from_lock_file()
outs = hub.execute_commands(commands)

# write outputs for debug
open("run.log", "w").write("\n".join(outs))
