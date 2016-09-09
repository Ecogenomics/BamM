from distutils.core import setup
import sys
from subprocess import check_call
from os.path import join, abspath, exists
from os import chdir, getcwd, rename, remove

xtra_opts = {"--with-libcfu-inc":"libcfu headers at this location",
             "--with-libcfu-lib":"libcfu library at this location",
             "--with-libhts-inc":"htslib headers at this location",
             "--with-libhts-lib":"htslib library at this location"}

if '-help' not in sys.argv and \
   '-h' not in sys.argv and \
   '--help' not in sys.argv and \
   '--h' not in sys.argv:
    if 'sdist' not in sys.argv:
        # grab extra configuration arguments
        configure_args = []
        for opt in xtra_opts.keys():
            try:
                opt_idx = sys.argv.index(opt)
                configure_args.append(opt+"="+abspath(sys.argv[opt_idx+1]))
                del sys.argv[opt_idx+1]
                sys.argv.remove(opt)
            except ValueError:
                pass

        # configure and make the c portion of the program
        cur_dir = getcwd()
        chdir(join(cur_dir, 'c'))
        check_call([join(getcwd(), "configure")] + configure_args)
        check_call(['make','clean'])
        check_call(['make'])
        chdir(cur_dir)
        # move the compiled library into the bamm folder so it's
        # picked up by distutils
        rename(join('c', 'libBamM.a'), join('bamm', 'libBamM.a'))

else:
    print
    print "Embedded C options (for building libPMBam.a) USE: --OPTION<space>PATH"
    for opt in xtra_opts.keys():
        print "  %s  %s"%(opt,xtra_opts[opt])
    print

exec(open('bamm/version.py').read()) # loads __version__

# return to regular viewing
setup(
    name='BamM',
    version=__version__,
    author='Michael Imelfort',
    author_email='mike@mikeimelfort.com',
    packages=['bamm'],
    scripts=['bin/bamm', 'bin/bamFlags'],
    license='LGPLv3',
    description='BamM - working with the BAM',
    long_description=open('README.md').read(),
    install_requires=["numpy >= 1.6.1"],
    package_data={'bamm' : ['libBamM.a']},
)

# remove the library
#to_remove = join('bamm', 'libBamM.a') 
#if exists(to_remove):
#    remove(to_remove)
