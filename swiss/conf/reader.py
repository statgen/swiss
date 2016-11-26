import yaml

def read_conf(fpath):
  with open(fpath,"rt") as fp:
    return yaml.safe_load(fp)

