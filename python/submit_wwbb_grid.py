#!/bin/env/python

from argparse import ArgumentParser
import glob
import re
import subprocess
import sys

def main() :
    parser = ArgumentParser(description = 'runWWbb grid submission')
    add_arg = parser.add_argument
    add_arg('-f', '--input-files', required = True, nargs = '*', help = 'input file with datasets, can specify more than one')
    add_arg('-v', '--verbose', action = 'store_true', default = False)
    add_arg('-p', '--pattern', default = '.*', help = 'grep pattern to select datasets')
    add_arg('-n', '--nickname', required = True, help = 'grid nickname, for naming output DS')
    add_arg('-t', '--tag', required = True, help = 'grid tag to add to dataset')
    add_arg('--destSE', default = 'SLACXRD_SCRATCHDISK', help = 'replicate output dataset to specified site')
    add_arg('--lumi', default = 35, help = 'set lumi (in fb) for job')
    add_arg('--suffix', default = '', help = 'set a job suffix')
    add_arg('--hh', action = 'store_true', default = False, help = 'set hh signal flag')
    add_arg('--bjet-eff', action = 'store_true', default = False, help = 'use bjet eff emulation')
    add_arg('--nFilesPerJob', default = '5',  help = 'prun option')
    add_arg('--nGBPerJob', help = 'prun option')
    add_arg('--group-role', action = 'store_true', help = 'submit jobs with group role')
    args = parser.parse_args()

    if args.nFilesPerJob and args.nGBPerJob :
        print parser.error("prun does not allow to specify both options '--nFilesPerJob' and '--nGBPerJob' at the same time")
    input_files = args.input_files

    print "Submitting {}\ninput file: {}\npattern:    {}".format(args.tag, input_files, args.pattern)
    for input_file in input_files :
        with open(input_file) as lines :
            lines = [l.strip() for l in lines if is_interesting_line(line=l, regexp=args.pattern)]
            for line in lines :
                inDS = line

                mc_type = "mc15_13TeV."
                if "mc16" in inDS :
                    mc_type = "mc16_13TeV."
                dsid = inDS[inDS.find(mc_type) + len(mc_type) : inDS.find(mc_type) + len(mc_type) + 6]
                sample = determine_sample_name(inDS)
                out_ds_suffix = 'nt'
                outDS = determine_outdataset_name(input_dataset_name = inDS, nt_tag = args.tag,
                            use_group = args.group_role, nickname = args.nickname,
                            prun_suffix = '_' + out_ds_suffix)


                # grid command
                grid_command = "./bash/gridScript.sh %IN "
                grid_command += (' --lumi ' + str(args.lumi))
                #grid_command += (' --suffix ' + args.suffix)
                grid_command += ('' if not args.hh else ' --hh')
                grid_command += ('' if not args.bjet_eff else ' --bjet-eff')
                grid_command += (' --dsid %s' % dsid)
                grid_command += (' --skip-maps') # skip use of local xsec and sumw map files for grid submission

                print "WARNING FORCING 5000 EVENTS TO RUN"
                print "WARNING FORCING 5000 EVENTS TO RUN"
                print "WARNING FORCING 5000 EVENTS TO RUN"
                print "WARNING FORCING 5000 EVENTS TO RUN"
                print "WARNING FORCING 5000 EVENTS TO RUN"
                grid_command += (' -n 5000')

                line_break = ('_'*90)
                print "\n{}\nInput {}\nOutput {}\nCommand {}".format(line_break, inDS, outDS, grid_command)

                exp_output_filename = "wwbb_truth_%s.root" % (dsid) #, args.suffix)

                # prun command
                prun_command = ('prun --exec "' + grid_command + '" --useRootCore --tmpDir /tmp ')
                prun_command += (' --inDS {} --outDS {}'.format(inDS, outDS))
                prun_command += (' --inTarBall=area.tgz --extFile "*.so,*.root" --match "*root*"')
                prun_command += (' --safetySize=600')
                prun_command += (' --outputs "{0}:{1}"'.format(out_ds_suffix, exp_output_filename))
                prun_command += (' --nFilesPerJob={}'.format(args.nFilesPerJob) if args.nFilesPerJob else '5')
                prun_command += (' --nGBPerJob={}'.format(args.nGBPerJob) if args.nGBPerJob else '')
                prun_command += (' --rootVer=6.04/16 --cmtConfig=x86_64-slc6-gcc49-opt')
                prun_command += ('' if not args.group_role else ' --official --voms atlas:/atlas/phys-susy/Role=production')
                #prun_command += (' --destSE=' + (args.destSE if not args.group_role else ','.join(args.destSE])))

                # execute prun
                if args.verbose : print prun_command
                subprocess.call(prun_command, shell = True)
                

def determine_outdataset_name(input_dataset_name, nt_tag, use_group, nickname, prun_suffix = "susyNt.root") :
    prefix = 'group.phys-susy.' if use_group else "user.%s." % nickname
    output_ds_name = prefix + re.sub('/', '', input_dataset_name) + '_' + nt_tag + '/'
    for i in range(5) :
        output_ds_name = re.sub('DAOD_TRUTH%d' % i, 'WWBB_TRUTH', output_ds_name)
    output_ds_name = re.sub('merge\.', '', output_ds_name)
    if output_ds_name.count('group.phys-susy.') > 1 :
        output_ds_name = output_ds_name.replace('group.phys-susy.', '', 1)
    max_ds_len = 132
    if len(output_ds_name + prun_suffix + '/') > max_ds_len :
        tags_to_keep = "_.*_%s" % nt_tag
        regex = "\.WWBB_TRUTH\.(?P<other_tags>.*)%s" % tags_to_keep
        match = re.search(regex, output_ds_name)
        if match :
            output_ds_name = output_ds_name.replace(match.group('other_tags'), '')
    output_ds_name = output_ds_name.replace('__', '_').replace('..', '.').replace('._', '.')
    return output_ds_name

def determine_sample_name(input_dataset_name = '') :
    sample = re.sub('merge.*', '', input_dataset_name)
    sample = re.sub('mc15_13TeV\.[0-9]*\.', '', sample)
    return sample

def is_interesting_line(line='', regexp='') :
    "interestin line = non-comment, non-empty, one name, matching selection"
    line = line.strip()
    tokens = line.split()
    return (len(line)>0 and not line.startswith('#') and len(tokens)==1 and re.search(regexp, line))


if __name__ == '__main__' :
    main()
