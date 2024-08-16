import subprocess
import sys
import os


def clean_pip():
    subprocess.run(
        [sys.executable, '-m', 'pip', 'uninstall', '-y', 'openseespy'])
    subprocess.run(
        ['rm', '-fr', 'build', 'dist', 'openseespy.egg-info'])


def getSetupCode(version, macversion):
    return f"""
import setuptools

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
    name="openseespy",
    version="{version}",
    author="Minjie Zhu",
    author_email="zhum@oregonstate.edu",
    description="A OpenSeesPy package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zhuminjie/openseespy",
    packages=setuptools.find_packages(),
    package_data={{
        '': ['LICENSE.md',
             '*.dat',
             '*.at2', ],
    }},
    license='LICENSE.md',
    classifiers=[
        "Programming Language :: Python",
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows',
        "Operating System :: MacOS :: MacOS X"],
    platforms=[
        "Linux",
        'Windows',
        'Mac'],
    install_requires=[
        f'openseespywin>={version}; platform_system=="Windows"',
        f'openseespylinux>={version}; platform_system=="Linux"',
        f'openseespymac>={macversion}; platform_system=="Darwin"',
    ],
    python_requires='>=3.6',
    zip_safe=False)
"""


def build_pip(version, macversion, use_zip):

    # change to script's folder
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    # clean folders
    clean_pip()

    # update tools
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--upgrade',
                    'setuptools', 'wheel', 'twine'],
                   check=True)

    # create setup.py
    with open('setup.py', 'w') as fd:
        fd.write(getSetupCode(version, macversion))

    # compile wheel
    if use_zip:
        subprocess.run(
            [sys.executable, 'setup.py',
                'bdist', '--format=zip'],
            check=True)
    else:
        subprocess.run([sys.executable, 'setup.py', 'bdist_wheel'],
                       check=True)


def upload_pip():
    # upload
    subprocess.run(
        [sys.executable, '-m', 'twine', 'upload', '--repository',
         'openseespy', 'dist/*'],
        check=True)


def upload_pip_test():
    # upload
    subprocess.run([sys.executable, '-m', 'twine', 'upload',
                    '--repository', 'testopenseespy', 'dist/*'],
                   check=True)


# python3 build_pip.py build version macversion
# python3 build_pip.py upload-test
# python3 build_pip.py upload
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print('build_pip cmd')
        exit()

    if sys.argv[1] == 'build':
        if len(sys.argv) < 4:
            print(
                'buld_pip build version macversion')
            exit()

        use_zip = False
        version = sys.argv[2]
        macversion = sys.argv[3]

        build_pip(version, macversion, use_zip)

    elif sys.argv[1] == 'upload-test':
        upload_pip_test()

    elif sys.argv[1] == 'upload':
        upload_pip()
