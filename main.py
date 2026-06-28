from importlib.metadata import version, PackageNotFoundError

def define_env(env):
    def get_version():
        try:
            return version("viewer3d")
        except PackageNotFoundError:
            return "dev"

    env.variables["version"] = get_version()
