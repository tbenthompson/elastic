from distutils.core import setup

def setup_package():
    setup(
        name = "elastic",
        version = "0.1",
        description = "A solver for the elastic boundary integral equations",
        author = "T. Ben Thompson",
        author_email = "t.ben.thompson@gmail.com",
        packages = ["elastic"]
    )

if __name__ == "__main__":
    setup_package()
