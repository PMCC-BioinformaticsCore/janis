#
# Utils to help out yaml dump
#


# Provide better formatting of multiline literals
# from: https://github.com/cloudlinux/kuberdock-platform/blob/master/kubedock/kapi/apps.py
def str_presenter(dumper, data):
    # check for multiline strings
    if len(data.splitlines()) == 1 and data[-1] == "\n":
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style=">")
    if len(data.splitlines()) > 1:
        return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
    return dumper.represent_scalar("tag:yaml.org,2002:str", data.strip())
