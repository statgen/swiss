import yaml

def write_conf(obj,fpath):
  with open(fpath,"wt") as fp:
    yaml.dump(obj,fp,default_flow_style=False,indent=2)

