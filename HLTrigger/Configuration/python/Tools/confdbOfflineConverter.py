#! /usr/bin/env python
import sys, os
import os.path
import tempfile
import urllib.request
import shutil
import subprocess
import atexit
from collections import Counter

class OfflineConverter:

    # the machine aliases and interfaces for the *online* database are
    #   cmsonr1-s.cms, cmsonr2-s.cms, cmsonr3-s.cms
    #   cmsonr1-v.cms, cmsonr2-v.cms, cmsonr3-v.cms
    # but the -s and -v interfaces resolve to the same hosts.
    # The actual machines and interfaces are
    #   CMSRAC11-S.cms, CMSRAC12-S.cms, CMSRAC21-S.cms
    #   CMSRAC11-V.cms, CMSRAC12-V.cms, CMSRAC21-V.cms

    # the possible machines and interfaces for the *offline* database are
    #   cmsr1-s.cms, cmsr2-s.cms, cmsr3-s.cms
    #   cmsr1-v.cms, cmsr2-v.cms, cmsr3-v.cms
    # but the -s and -v interfaces resolve to the same hosts
    # The actual machines and interfaces are
    #   itrac50011-s.cern.ch, itrac50063-s.cern.ch, itrac50078-s.cern.ch
    #   itrac50011-v.cern.ch, itrac50063-v.cern.ch, itrac50078-v.cern.ch

    databases = {}
    databases['v1'] = {}
    databases['v1']['offline'] = ( '-t', 'oracle', '-h', 'cmsr1-s.cern.ch',        '-d', 'cms_cond.cern.ch',      '-u', 'cms_hltdev_reader', '-s', 'convertMe!' )
    databases['v1']['hltdev']  = databases['v1']['offline']     # for backwards compatibility
    databases['v1']['online']  = ( '-t', 'oracle', '-h', 'cmsonr1-s.cms',          '-d', 'cms_rcms.cern.ch',      '-u', 'cms_hlt_r',         '-s', 'convertMe!' )
    databases['v1']['adg']     = ( '-t', 'oracle', '-h', 'cmsr1-s.cern.ch',        '-d', 'cms_cond.cern.ch',      '-u', 'cms_hlt_gui_r',     '-s', 'convertMe!' )
    databases['v1']['orcoff']  = databases['v1']['adg']         # for backwards compatibility
    databases['v3'] = {}
    databases['v3']['run2'] = ( '-t', 'oracle', '-h', 'cmsr1-s.cern.ch,cmsr2-s.cern.ch,cmsr3-s.cern.ch',        '-d', 'cms_hlt.cern.ch',      '-u', 'cms_hlt_gdr_r',     '-s', 'convertMe!' )
    databases['v3']['run3'] = ( '-t', 'oracle', '-h', 'cmsr1-s.cern.ch,cmsr2-s.cern.ch,cmsr3-s.cern.ch',        '-d', 'cms_hlt.cern.ch',      '-u', 'cms_hlt_v3_r',     '-s', 'convertMe!' )
    databases['v3']['dev'] = ( '-t', 'oracle', '-h', 'cmsr1-s.cern.ch,cmsr2-s.cern.ch,cmsr3-s.cern.ch',        '-d', 'cms_hlt.cern.ch',      '-u', 'cms_hlt_gdrdev_r',     '-s', 'convertMe1!' )
    databases['v3']['online']  = ( '-t', 'oracle', '-h', 'cmsonr1-s.cms',          '-d', 'cms_rcms.cern.ch',      '-u', 'cms_hlt_gdr_r',     '-s', 'convertMe!' )
    databases['v3']['adg']     = ( '-t', 'oracle', '-h', 'cmsonr1-adg1-s.cern.ch', '-d', 'cms_orcon_adg.cern.ch', '-u', 'cms_hlt_gdr_r',     '-s', 'convertMe!' )
    databases['v3-beta'] = dict(databases['v3'])
    databases['v3-test'] = dict(databases['v3'])
    databases['v2'] = dict(databases['v3'])
    #old converter can only handle a single host so we modify the params accordingly
    for dbkey in databases['v2']:
        dbparams  = databases['v2'][dbkey]
        if dbparams[3]=='cmsr1-s.cern.ch,cmsr2-s.cern.ch,cmsr3-s.cern.ch':
            databases['v2'][dbkey] = dbparams[0:3]+('cmsr1-s.cern.ch',)+dbparams[4:]

    @staticmethod
    def CheckTempDirectory(dir):
        dir = os.path.realpath(dir)
        if not os.path.isdir(dir):
            try:
                os.makedirs(dir)
            except:
                return None
        return dir


    def __init__(self, version = 'v3', database = 'run3', url = None, verbose = False):
        self.verbose = verbose
        self.version = version
        self.baseDir = '/afs/cern.ch/user/c/confdb/www/%s/lib' % version
        self.baseUrl = 'https://confdb.web.cern.ch/confdb/%s/lib' % version
        self.jars    = ( 'ojdbc8.jar', 'cmssw-evf-confdb-converter.jar' )
        if version=='v2':
            #legacy driver for run2 gui
            self.jars = ( 'ojdbc6.jar', 'cmssw-evf-confdb-converter.jar' )
        self.workDir = ''

        # check the schema version
        if version not in self.databases:
            # unsupported database version
            sys.stderr.write( "ERROR: unsupported database version \"%s\"\n" % version)

        # check the database
        if database in self.databases[version]:
            # load the connection parameters for the given database
            self.connect = self.databases[version][database]
        else:
            # unsupported database
            sys.stderr.write( "ERROR: unknown database \"%s\" for version \"%s\"\n" % (database, version))
            sys.exit(1)

        # check for a custom base URL
        if url is not None:
            self.baseUrl = url

        # try to read the .jar files from AFS, or download them
        if os.path.isdir(self.baseDir) and all(os.path.isfile(self.baseDir + '/' + jar) for jar in self.jars):
            # read the .jar fles from AFS
            self.workDir = self.baseDir
        else:
            # try to use $CMSSW_BASE/tmp
            self.workDir = OfflineConverter.CheckTempDirectory(os.environ['CMSSW_BASE'] + '/tmp/confdb')
            if not self.workDir:
                # try to use $TMP
                self.workDir = OfflineConverter.CheckTempDirectory(os.environ['TMP'] + '/confdb')
            if not self.workDir:
                # create a new temporary directory, and install a cleanup callback
                self.workDir = tempfile.mkdtemp()
                atexit.register(shutil.rmtree, self.workDir)
            # download the .jar files
            for jar in self.jars:
                # check if the file is already present
                if os.path.exists(self.workDir + '/' + jar):
                    continue
                # download to a temporay name and use an atomic rename (in case an other istance is downloading the same file
                handle, temp = tempfile.mkstemp(dir = self.workDir, prefix = jar + '.')
                os.close(handle)
                urllib.request.urlretrieve(self.baseUrl + '/' + jar, temp)
                if not os.path.exists(self.workDir + '/' + jar):
                    os.rename(temp, self.workDir + '/' + jar)
                else:
                    os.unlink(temp)

        # setup the java command line and CLASSPATH
        if self.verbose:
            sys.stderr.write("workDir = %s\n" % self.workDir)
        # use non-blocking random number source /dev/urandom (instead of /dev/random), see:
        #   http://blockdump.blogspot.fr/2012/07/connection-problems-inbound-connection.html
        # deal with timezone region not found
        #   http://stackoverflow.com/questions/9156379/ora-01882-timezone-region-not-found
        # increase the thread stack size from the default of 1 MB to work around java.lang.StackOverflowError errors, see
        #   man java
        self.javaCmd = ( 'java', '-cp', ':'.join(self.workDir + '/' + jar for jar in self.jars), '-Djava.security.egd=file:///dev/urandom', '-Doracle.jdbc.timezoneAsRegion=false', '-Xss32M', 'confdb.converter.BrowserConverter' )


    def query(self, *args):
        args = self.javaCmd + self.connect + args
        if self.verbose:
            sys.stderr.write("\n" + ' '.join(args) + "\n\n" )
        sub = subprocess.Popen(
            args,
            stdin = None,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            shell = False,
            universal_newlines = True )
        return sub.communicate()

def help():
    sys.stdout.write("""Usage: %s OPTIONS

        --v1|--v2|--v3|--v3-beta|--v3-test  (specify the ConfDB version [default: v3])

        --run3|--run2|--dev|--online|--adg    (specify the target db [default: run3], online will only work inside p5 network)

        Note that for v1
            --orcoff  is a synonim of --adg
            --offline is a synonim of --hltdev

        --configId <id>             (specify the configuration by id)
        --configName <name>         (specify the configuration by name)
        --runNumber <run>           (specify the configuration by run number)
          [exactly one of --configId OR --configName OR --runNumber is required]

        --cff                       (retrieve configuration *fragment*)
        --input <f1.root[,f2.root]> (insert PoolSource with specified fileNames)
        --input <files.list>        (read a text file which lists input ROOT files)
        --output <out.root>         (insert PoolOutputModule w/ specified fileName)
        --nopsets                   (exclude all globale psets)
        --noedsources               (exclude all edsources)
        --noes                      (exclude all essources *and* esmodules)
        --noessources               (exclude all essources)
        --noesmodules               (exclude all esmodules)
        --noservices                (exclude all services)
        --nooutput                  (exclude all output modules)
        --nopaths                   (exclude all paths [+=referenced seqs&mods])
        --nosequences               (don't define sequences [+=referenced s&m])
        --nomodules                 (don't define modules)
        --psets <pset1[,pset2]>     (include only specified global psets)
        --psets <-pset1[,-pset2]>   (include all global psets but the specified)
        --essources <ess1[,ess2]>   (include only specified essources)
        --essources <-ess1[,-ess2]> (include all essources but the specified)
        --esmodules <esm1[,esm2]>   (include only specified esmodules)
        --esmodules <-esm1[,-esm2]> (include all esmodules but the specified)
        --services <svc1[,svc2]>    (include only specified services)
        --services <-svc1[,-svc2]>  (include all services but the specified)
        --paths <p1[,p2]>           (include only specified paths)
        --paths <-p1[,-p2]>         (include all paths but the specified)
        --streams <s1[,s2]>         (include only specified streams)
        --datasets <d1[,d2]>        (include only specified datasets)
        --sequences <s1[,s2]>       (include sequences, referenced or not!)
        --modules <p1[,p2]>         (include modules, referenced or not!)
        --blocks <m1::p1[,p2][,m2]> (generate parameter blocks)

        --verbose                   (print additional details)
""")


def main():
    args = sys.argv[1:]
    version = 'v3'
    db      = 'run3'
    verbose = False

    if not args:
        help()
        sys.exit(1)

    if '--help' in args or '-h' in args:
        help()
        sys.exit(0)

    if '--verbose' in args:
        verbose = True
        args.remove('--verbose')

    arg_count = Counter(args)
    db_count = arg_count['--v1'] + arg_count['--v2'] + arg_count['--v3'] + arg_count['--v3-beta'] + arg_count['--v3-test']
    if db_count>1:
        sys.stderr.write( 'ERROR: conflicting database version specifications: "--v1", "--v2", "--v3", "--v3-beta", and "--v3-test" are mutually exclusive options' )
        sys.exit(1)

    if '--v1' in args:
        version = 'v1'
        db      = 'offline'
        args.remove('--v1')

    if '--v2' in args:
        version = 'v2'
        db      = 'run2'
        args.remove('--v2')

    if '--v3' in args:
        version = 'v3'
        db      = 'run3'
        args.remove('--v3')

    if '--v3-beta' in args:
        version = 'v3-beta'
        db      = 'run3'
        args.remove('--v3-beta')

    if '--v3-test' in args:
        version = 'v3-test'
        db      = 'dev'
        args.remove('--v3-test')

    _dbs = {}
    _dbs['v1'] = [ '--%s' % _db for _db in OfflineConverter.databases['v1'] ] + [ '--runNumber' ]
    _dbs['v2'] = [ '--%s' % _db for _db in OfflineConverter.databases['v2'] ] + [ '--runNumber' ]
    _dbs['v3'] = [ '--%s' % _db for _db in OfflineConverter.databases['v3'] ] + [ '--runNumber'] 
    _dbs['v3-beta'] = [ '--%s' % _db for _db in OfflineConverter.databases['v3-beta'] ] + [ '--runNumber' ]
    _dbs['v3-test'] = [ '--%s' % _db for _db in OfflineConverter.databases['v3-test'] ] + [ '--runNumber' ]
    _dbargs = set(args) & set(sum(_dbs.values(), []))

    if _dbargs:
        if len(_dbargs) > 1:
            sys.stderr.write( "ERROR: too many database specifications: \"" + "\", \"".join( _dbargs) + "\"\n" )
            sys.exit(1)

        _arg = _dbargs.pop()
        db   = _arg[2:]
        if db == 'runNumber':
            db = 'adg'
        else:
            args.remove(_arg)

        if not db in OfflineConverter.databases[version]:
            sys.stderr.write( "ERROR: database version \"%s\" incompatible with specification \"%s\"\n" % (version, db) )
            sys.exit(1)

    converter = OfflineConverter(version = version, database = db, verbose = verbose)
    out, err = converter.query( * args )
    if 'ERROR' in err:
        sys.stderr.write( "%s: error while retriving the HLT menu\n\n%s\n\n" % (sys.argv[0], err) )
        sys.exit(1)
    else:
        sys.stdout.write( out )


if __name__ == "__main__":
    main()
