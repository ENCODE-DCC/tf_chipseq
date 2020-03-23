#!/usr/bin/env python2
'''Instantiate the ENCODE ChIP-seq workflow'''

import sys
import logging
import re
import dxpy
import time
import pprint

EPILOG = '''Notes:

Examples:
    # Build blank TF workflow from fastq to peaks
    %(prog)s --target tf --name "ENCODE TF ChIP-seq (no reference)" --outf "/ChIP-seq/"

    # Build blank histone workflow from fastq to peaks
    %(prog)s --target histone --name "ENCODE Histone ChIP-seq (no reference)" --outf "/ChIP-seq/"

    # Build a pre-configured GRCh38 histone workflow, requiring only data to run
    %(prog)s --target histone \\
    --name "ENCODE Histone ChIP-seq (GRCh38)" \\
    --chrom_sizes "ENCODE Reference Files:/GRCh38/GRCh38_EBV.chrom.sizes" \\
    --genomesize hs \\
    --reference "ENCODE Reference Files:/GRCh38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fa.gz" \\
    --outf "/ChIP-seq/"

    # Build and run a complete hg19 TF workflow, specifying all inputs.
    %(prog)s --target tf \\
    --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \\
    --genomesize hs \\
    --reference "ENCODE Reference Files:/hg19/male.hg19.tar.gz" \\
    --blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \\
    --outf "ENCSR464DKE-hCTCF-chr21" \\
    --title "ENCSR464DKE-hCTCF-chr21" \\
    --rep1 "/ChIP-seq/test_data/ENCSR464DKE-hCTCF/R1-ENCFF921SED.chr21.fq.gz" \\
    --rep2 "/ChIP-seq/test_data/ENCSR464DKE-hCTCF/R2-ENCFF812KOM.chr21.fq.gz" \\
    --ctl1 "/ChIP-seq/test_data/ENCSR464DKE-hCTCF/C1-ENCFF690VPV.chr21.fq.gz" \\
    --ctl2 "/ChIP-seq/test_data/ENCSR464DKE-hCTCF/C2-ENCFF357TLV.chr21.fq.gz" \\
    --yes

    # Build and run a complete hg19 TF workflow, with a unary control.
    %(prog)s --target tf \\
    --chrom_sizes "ENCODE Reference Files:/hg19/male.hg19.chrom.sizes" \\
    --genomesize hs \\
    --reference "ENCODE Reference Files:/hg19/male.hg19.tar.gz" \\
    --blacklist "ENCODE Reference Files:/hg19/blacklists/wgEncodeDacMapabilityConsensusExcludable.bed.gz" \\
    --outf "ENCSR000EEB-hMAFK-chr21" \\
    --title "ENCSR000EEB-hMAFK-chr21" \\
    --rep1 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/R1-ENCFF000XTT.chr21.fq.gz" \\
    --rep2 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/R2-ENCFF000XTU.chr21.fq.gz" \\
    --ctl1 "/ChIP-seq/test_data/ENCSR000EEB-hMAFK/C1-ENCFF000XSJ.chr21.fq.gz" \\
    --yes

    # Build and run a complete mm10 histone workflow, specifying all inputs.
    %(prog)s --target histone \\
    --chrom_sizes "ENCODE Reference Files:/mm10/male.mm10.chrom.sizes" \\
    --genomesize mm \\
    --reference "ENCODE Reference Files:/mm10/male.mm10.tar.gz" \\
    --outf "ENCSR087PLZ-mH3K9ac-chr19" \\
    --title "ENCSR087PLZ-mH3K9ac-chr19" \\
    --rep1 "/ChIP-seq/test_data/ENCSR087PLZ-mH3K9ac/R1-ENCFF560GLI.chr19.fq.gz" \\
    --rep2 "/ChIP-seq/test_data/ENCSR087PLZ-mH3K9ac/R2-ENCFF891NNX.chr19.fq.gz" \\
    --ctl1 "/ChIP-seq/test_data/ENCSR087PLZ-mH3K9ac/C1-ENCFF069WCH.chr19.fq.gz" \\
    --ctl2 "/ChIP-seq/test_data/ENCSR087PLZ-mH3K9ac/C2-ENCFF101KOM.chr19.fq.gz" \\
    --yes
'''

WF = {
    'default': {
        'wf_name': 'chip_seq',
        'wf_title': 'ChIP-seq',
        'wf_description': 'ENCODE ChIP-seq Analysis Pipeline',
        'run_idr': True
    },
    'histone': {
        'wf_name': 'histone_chip_seq',
        'wf_title': 'Histone ChIP-seq',
        'wf_description': 'ENCODE histone ChIP-seq Analysis Pipeline',
        'run_idr': False
    },
    'tf': {
        'wf_name': 'tf_chip_seq',
        'wf_title': 'TF ChIP-seq',
        'wf_description': 'ENCODE TF ChIP-seq Analysis Pipeline',
        'run_idr': True
    }
}

DEFAULT_APPLET_PROJECT = dxpy.WORKSPACE_ID
DEFAULT_OUTPUT_PROJECT = dxpy.WORKSPACE_ID
DEFAULT_OUTPUT_FOLDER = '/analysis_run'

MAPPING_APPLET_NAME = 'encode_map'
FILTER_QC_APPLET_NAME = 'filter_qc'
XCOR_APPLET_NAME = 'xcor'
XCOR_ONLY_APPLET_NAME = 'xcor_only'
SPP_APPLET_NAME = 'spp'
POOL_APPLET_NAME = 'pool'
PSEUDOREPLICATOR_APPLET_NAME = 'pseudoreplicator'
ENCODE_SPP_APPLET_NAME = 'encode_spp'
ENCODE_MACS2_APPLET_NAME = 'encode_macs2'
# IDR_APPLET_NAME='idr'
IDR2_APPLET_NAME = 'idr2'
ENCODE_IDR_APPLET_NAME = 'encode_idr'
OVERLAP_PEAKS_APPLET_NAME = 'overlap_peaks'
ACCESSION_ANALYSIS_APPLET_NAME = 'accession_analysis'

APPLETS = {}


def get_args():
    import argparse
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    def t_or_f(arg):
        ua = str(arg).upper()
        if ua == 'TRUE'[:len(ua)]:
            return True
        elif ua == 'FALSE'[:len(ua)]:
            return False
        else:
            assert not (True or False), "Cannot parse %s to boolean" % (arg)

    parser.add_argument(
        '--target',
        help="ChIP target type (histone or tf)",
        required=True)
    parser.add_argument(
        '--debug',
        help="Print debug messages and hold jobs for ssh",
        default=False, action='store_true')
    parser.add_argument(
        '--reference',
        help="Reference tar to map to")
    parser.add_argument(
        '--chrom_sizes',
        help="chrom.sizes file for bedToBigBed")
    parser.add_argument(
        '--genomesize',
        help="Genome size string for MACS2, e.g. mm or hs")
    parser.add_argument(
        '--narrowpeak_as',
        help=".as file for bed to bigbed",
        default='ENCODE Reference Files:narrowPeak.as')
    parser.add_argument(
        '--gappedpeak_as',
        help=".as file for bed to bigbed",
        default='ENCODE Reference Files:gappedPeak.as')
    parser.add_argument(
        '--broadpeak_as',
        help=".as file for bed to bigbed",
        default='ENCODE Reference Files:broadPeak.as')
    parser.add_argument(
        '--rep1',
        help="Replicate 1 fastq or tagAlign",
        default=None, nargs='*')
    parser.add_argument(
        '--rep2',
        help="Replicate 2 fastq or tagAlign",
        default=None, nargs='*')
    parser.add_argument(
        '--ctl1',
        help="Control for replicate 1 fastq or tagAlign",
        default=None, nargs='*')
    parser.add_argument(
        '--ctl2',
        help="Control for replicate 2 fastq or tagAlign",
        default=None,  nargs='*')
    parser.add_argument(
        '--unary_control',
        help="Force one control for both reps",
        default=False, action='store_true')
    parser.add_argument(
        '--simplicate_experiment',
        help="Force single replicate (rep1)",
        default=False, action='store_true')
    parser.add_argument(
        '--outp',
        help="Output project name or ID",
        default=DEFAULT_OUTPUT_PROJECT)
    parser.add_argument(
        '--outf',
        help="Output folder name or ID",
        default=DEFAULT_OUTPUT_FOLDER)
    parser.add_argument(
        '--use_existing_folders',
        help="Reuse existing folders even if results have already been saved there",
        default=False, action='store_true')
    parser.add_argument(
        '--name',
        help="Name for new workflow")
    parser.add_argument(
        '--title',
        help="Title for new workflow")
    parser.add_argument(
        '--description',
        help="Description for new workflow")
    parser.add_argument(
        '--applets',
        help="Name of project containing applets",
        default=DEFAULT_APPLET_PROJECT)
    parser.add_argument(
        '--nomap',
        help='Given tagAligns, skip to peak calling',
        default=False, action='store_true')
    parser.add_argument(
        '--maponly',
        help='Given fastqs, only map and calculate xcor but no peaks',
        default=False, action='store_true')
    parser.add_argument(
        '--rep1pe',
        help='Specify if rep1 is PE (required only if --nomap)',
        type=t_or_f, default=None)
    parser.add_argument(
        '--rep2pe',
        help='Specify if rep2 is PE (required only if --nomap)',
        type=t_or_f, default=None)
    parser.add_argument(
        '--blacklist',
        help="Blacklist to filter IDR peaks")
    parser.add_argument(
        '--yes',
        help='Run the workflow',
        default=False, action='store_true')
    parser.add_argument(
        '--spp_version',
        help="Version string for spp",
        default="1.14")
    parser.add_argument(
        '--spp_instance',
        help="Override the spp applet instance type",
        default=None)
    parser.add_argument(
        '--pipeline_version',
        help="Version string for ENCODE pipeline",
        default="1.2")
    parser.add_argument(
        '--accession',
        help='Automatically accession the results to the ENCODE Portal',
        default=False, action='store_true')
    parser.add_argument('--fqcheck', help="If --accession, check that analysis is based on latest fastqs on ENCODEd", type=t_or_f, default=None)
    parser.add_argument('--skip_control', help="If --accession, accession no control files or metadata", type=t_or_f, default=None)
    parser.add_argument('--force_patch', help="Force patching metadata for existing files", type=t_or_f, default=None)
    parser.add_argument('--scrub', help="Also produce bams scrubbed of sequence information", type=t_or_f, default=None)
    parser.add_argument('--fragment_length',
                        type=int,
                        help="Instead of calculating fragment length from xcor, use this fragment length",
                        default=None)
    # parser.add_argument('--idr',     help='Report peaks with and without IDR analysis',                 default=False, action='store_true')
    # parser.add_argument('--idronly',  help='Only report IDR peaks', default=None, action='store_true')
    # parser.add_argument('--idrversion', help='Version of IDR to use (1 or 2)', default="2")

    args = parser.parse_args()

    global DEBUG
    DEBUG = args.debug
    if DEBUG:
        logging.basicConfig(
            format='%(levelname)s:%(message)s',
            level=logging.DEBUG)
        logging.debug("Debug logging ON")
    else:  # use the defaulf logging level
        logging.basicConfig(
            format='%(levelname)s:%(message)s')

    logging.debug("rep1 is: %s" % (args.rep1))

    return args


def blank_workflow(args):
    return


def map_and_filter(infile, args):
    if not infile:
        return {None}
    stages = {None}
    return stages


def call_peaks(expvsctl, args):
    if not expvsctl:
        return {None}
    stages = {None}
    return stages


def resolve_project(identifier, privs='r'):
    project = dxpy.find_one_project(
        name=identifier,
        level='VIEW',
        name_mode='exact',
        return_handler=True,
        zero_ok=True)
    if project is None:
        try:
            project = dxpy.get_handler(identifier)
        except:
            logging.error(
                'Could not find a unique project with name or id %s'
                % (identifier))
            raise ValueError(identifier)
    logging.debug(
        'Project %s access level is %s'
        % (project.name, project.describe()['level']))
    if privs == 'w' and project.describe()['level'] == 'VIEW':
        logging.error('Output project %s is read-only' % (identifier))
        raise ValueError(identifier)
    return project


def create_folder(project, folder_name):
    if not folder_name.startswith('/'):
        folder_name = '/' + folder_name
    try:
        project.new_folder(folder_name, parents=True)
    except:
        logging.error(
            "Cannot create folder %s in project %s"
            % (folder_name, project.name))
        return None
    else:
        logging.info(
            "New folder %s created in project %s"
            % (folder_name, project.name))
        return folder_name


def resolve_folder(project, identifier):
    if not identifier.startswith('/'):
        identifier = '/' + identifier
    try:
        project.list_folder(identifier)
    except:
        return None
    else:
        return identifier


def resolve_file(identifier):
    logging.debug("resolve_file: %s" % (identifier))

    if not identifier:
        return None

    m = re.match(r'''^([\w\-\ \.]+):([\w\-\ /\.]+)''', identifier)
    if m:
        project_identifier = m.group(1)
        file_identifier = m.group(2)
    else:
        logging.debug("Defaulting to the current project")
        project_identifier = dxpy.WORKSPACE_ID
        file_identifier = identifier

    project = resolve_project(project_identifier)
    logging.debug("Got project %s" % (project.name))
    logging.debug("Now looking for file %s" % (file_identifier))

    m = re.match(r'''(^[\w\-\ /\.]+)/([\w\-\ \.]+)''', file_identifier)
    if m:
        folder_name = m.group(1)
        if not folder_name.startswith('/'):
            folder_name = '/' + folder_name
        recurse = False
        file_name = m.group(2)
    else:
        folder_name = '/'
        recurse = True
        file_name = file_identifier

    logging.debug(
        "Looking for file %s in folder %s" % (file_name, folder_name))

    try:
        file_handler = dxpy.find_one_data_object(
            name=file_name,
            folder=folder_name,
            project=project.get_id(),
            recurse=recurse,
            more_ok=False,
            zero_ok=False,
            return_handler=True)
    except dxpy.DXSearchError:
        logging.debug(
            '%s not found in project %s folder %s.  Trying as file ID'
            % (file_name, project.get_id(), folder_name))
        file_handler = None
    except:
        raise

    if not file_handler:
        try:
            file_handler = dxpy.DXFile(dxid=identifier, mode='r')
        except dxpy.DXError:
            logging.debug('%s not found as a dxid' % (identifier))
            logging.warning('Could not find file %s.' % (identifier))
            file_handler = None
        except:
            raise

    if file_handler:
        logging.info(
            "Resolved file identifier %s to %s"
            % (identifier, file_handler.get_id()))
        return file_handler
    else:
        logging.warning("Failed to resolve file identifier %s" % (identifier))
        return None


def find_applet_by_name(applet_name, applets_project_id):
    '''Looks up an applet by name in the project that holds tools.
      From Joe Dale's code.'''
    cached = '*'
    if (applet_name, applets_project_id) not in APPLETS:
        found = dxpy.find_one_data_object(
            classname="applet",
            name=applet_name,
            project=applets_project_id,
            zero_ok=False,
            more_ok=False,
            return_handler=True)
        APPLETS[(applet_name, applets_project_id)] = found
        cached = ''

    logging.info(
        cached + "Resolved applet %s to %s"
        % (applet_name, APPLETS[(applet_name, applets_project_id)].get_id()))
    return APPLETS[(applet_name, applets_project_id)]


def main():
    args = get_args()

    blank_workflow = not (args.rep1 or args.rep2 or args.ctl1 or args.ctl2)

    if not blank_workflow:
        assert args.rep1, "Reads are required for rep1"
        assert args.ctl1, "Reads are required for ctl1"
        assert not args.nomap or args.rep1pe is not None, "With --nomap, endedness of rep1 must be specified witn --rep1pe"
        assert not args.nomap or (not args.rep2 or args.rep2pe is not None), "With --nomap, endedness of rep2 must be specified with --rep2pe"

    if not args.target:
        target_type = 'default'  # default
    else:
        target_type = args.target.lower()
    if target_type not in WF.keys():
        logging.error('Target type %s is not recognized')
        sys.exit(2)

    output_project = resolve_project(args.outp, 'w')
    logging.debug('Found output project %s' % (output_project.name))
    applet_project = resolve_project(args.applets, 'r')
    logging.debug('Found applet project %s' % (applet_project.name))    

    existing_folder = resolve_folder(output_project, args.outf)
    if not existing_folder:
        output_folder = create_folder(output_project, args.outf)
    elif args.use_existing_folders:
        output_folder = existing_folder
    else:
        assert (existing_folder and args.use_existing_folders), 'Output folder %s exists but --use_existing_folders is %s' % (existing_folder, args.use_existing_folders)

    logging.debug('Using output folder %s' % (output_folder))

    workflow = dxpy.new_dxworkflow(
        name=args.name or WF[target_type]['wf_name'],
        title=args.title or WF[target_type]['wf_title'],
        description=args.description or WF[target_type]['wf_description'],
        project=output_project.get_id(),
        folder=output_folder,
        properties={'pipeline_version': str(args.pipeline_version)})

    unary_control = args.unary_control or (not blank_workflow and args.ctl2 is None)
    simplicate_experiment = args.simplicate_experiment or (args.rep1 and not args.rep2)

    if not args.genomesize:
        genomesize = None
    else:
        genomesize = args.genomesize
    if not args.chrom_sizes:
        chrom_sizes = None
    else:
        chrom_sizes = dxpy.dxlink(resolve_file(args.chrom_sizes))

    if not args.blacklist:
        blacklist = None
    else:
        blacklist = dxpy.dxlink(resolve_file(args.blacklist))

    run_idr = WF[target_type]['run_idr']

    if not args.nomap:
        # a "superstage" is just a dict with a name, name(s) of input files,
        # and then names and id's of stages that process that input
        # each superstage here could be implemented as a stage in a more
        # abstract workflow.  That stage would then call the various applets
        # that are separate
        # stages here.
        mapping_superstages = [  # the order of this list is important in that
            {'name': 'Rep1', 'input_args': args.rep1}
        ]
        if not simplicate_experiment:
            mapping_superstages.append(
                {'name': 'Rep2', 'input_args': args.rep2})
        mapping_superstages.append(
            {'name': 'Ctl1', 'input_args': args.ctl1})
        if not unary_control and not simplicate_experiment:
            mapping_superstages.append(
                {'name': 'Ctl2', 'input_args': args.ctl2})

        mapping_applet = find_applet_by_name(
            MAPPING_APPLET_NAME, applet_project.get_id())
        # mapping_output_folder = resolve_folder(
        #     output_project, output_folder + '/' + mapping_applet.name)
        mapping_output_folder = mapping_applet.name
        reference_tar = resolve_file(args.reference)
        filter_qc_applet = find_applet_by_name(
            FILTER_QC_APPLET_NAME, applet_project.get_id())
        filter_qc_output_folder = mapping_output_folder
        xcor_applet = find_applet_by_name(
            XCOR_APPLET_NAME, applet_project.get_id())
        xcor_output_folder = mapping_output_folder

        # in the first pass create the mapping stage id's so we can use JBOR's
        # to link inputs
        for mapping_superstage in mapping_superstages:
            superstage_name = mapping_superstage.get('name')
            mapped_stage_id = workflow.add_stage(
                mapping_applet,
                name='Map %s' % (superstage_name),
                folder=mapping_output_folder
            )
            mapping_superstage.update({'map_stage_id': mapped_stage_id})

        # in the second pass populate the stage inputs and build other stages
        rep1_stage_id = next(ss.get('map_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1')
        for mapping_superstage in mapping_superstages:
            superstage_name = mapping_superstage.get('name')
            superstage_id = mapping_superstage.get('map_stage_id')

            if mapping_superstage.get('input_args') or blank_workflow:
                mapping_stage_input = {}
                if superstage_name != "Rep1":
                    mapping_stage_input.update(
                        {'reference_tar': dxpy.dxlink(
                            {'stage': rep1_stage_id,
                             'inputField': 'reference_tar'})})
                else:
                    if args.reference:
                        mapping_stage_input.update(
                            {'reference_tar': dxpy.dxlink(
                                reference_tar.get_id())})
                if not blank_workflow:
                    for arg_index, input_arg in enumerate(mapping_superstage['input_args']): #read pairs assumed be in order read1,read2
                        reads = dxpy.dxlink(resolve_file(input_arg).get_id())
                        mapping_stage_input.update({'reads%d' %(arg_index+1): reads})
                # this is now done in the first pass loop above
                # mapped_stage_id = workflow.add_stage(
                #     mapping_applet,
                #     name='Map %s' %(superstage_name),
                #     folder=mapping_output_folder,
                #     stage_input=mapping_stage_input
                # )
                # mapping_superstage.update({'map_stage_id': mapped_stage_id})
                workflow.update_stage(superstage_id, stage_input=mapping_stage_input)

                filter_qc_stage_input = {
                    'input_bam': dxpy.dxlink({'stage': superstage_id, 'outputField': 'mapped_reads'}),
                    'paired_end': dxpy.dxlink({'stage': superstage_id, 'outputField': 'paired_end'})
                }
                if args.scrub is not None:
                    filter_qc_stage_input.update({'scrub': args.scrub})
                filter_qc_stage_id = workflow.add_stage(
                    filter_qc_applet,
                    name='Filter_QC %s' %(superstage_name),
                    folder=filter_qc_output_folder,
                    stage_input=filter_qc_stage_input
                )
                mapping_superstage.update({'filter_qc_stage_id': filter_qc_stage_id})

                xcor_stage_id = workflow.add_stage(
                    xcor_applet,
                    name='Xcor %s' %(superstage_name),
                    folder=xcor_output_folder,
                    stage_input={
                        'input_bam': dxpy.dxlink({'stage': filter_qc_stage_id, 'outputField': 'filtered_bam'}),
                        'paired_end': dxpy.dxlink({'stage': filter_qc_stage_id, 'outputField': 'paired_end'}),
                        'spp_version': args.spp_version
                    }
                )
                mapping_superstage.update({'xcor_stage_id': xcor_stage_id})

        exp_rep1_ta = dxpy.dxlink(
                    {'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
                     'outputField': 'tagAlign_file'})
        exp_rep1_cc = dxpy.dxlink(
                    {'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
                     'outputField': 'CC_scores_file'})
        rep1_paired_end = dxpy.dxlink(
                        {'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep1'),
                         'outputField': 'paired_end'})
        if not simplicate_experiment:
            exp_rep2_ta = dxpy.dxlink(
                        {'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
                         'outputField': 'tagAlign_file'})
            exp_rep2_cc = dxpy.dxlink(
                        {'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
                         'outputField': 'CC_scores_file'})
            rep2_paired_end = dxpy.dxlink(
                            {'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Rep2'),
                             'outputField': 'paired_end'})
        else:
            exp_rep2_ta = None
            exp_rep2_cc = None
            rep2_paired_end = None

        ctl_rep1_ta = dxpy.dxlink(
                    {'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Ctl1'),
                     'outputField': 'tagAlign_file'})
        if not unary_control and not simplicate_experiment:
            ctl_rep2_ta = dxpy.dxlink(
                        {'stage': next(ss.get('xcor_stage_id') for ss in mapping_superstages if ss['name'] == 'Ctl2'),
                         'outputField': 'tagAlign_file'})
        else:
            ctl_rep2_ta = None

    else:  # skipped the mapping, so just bring in the inputs from arguments
        if not blank_workflow:
            exp_rep1_ta = dxpy.dxlink(resolve_file(args.rep1[0]).get_id())
            exp_rep1_ta_desc = dxpy.describe(exp_rep1_ta)
            exp_rep1_mapping_analysis_id = dxpy.describe(exp_rep1_ta_desc['createdBy']['job'])['analysis']
            exp_rep1_mapping_analysis = dxpy.describe(exp_rep1_mapping_analysis_id)
            rep1_xcor_stage_description = next(
                stage
                for stage in exp_rep1_mapping_analysis.get('stages')
                if stage['execution']['executableName'] == 'xcor')
            exp_rep1_cc = rep1_xcor_stage_description['execution']['output']['CC_scores_file']
            if args.rep1pe is None:
                print("Inferring rep1 PE-ness from analysis")
                rep1_paired_end = rep1_xcor_stage_description['execution']['output']['paired_end']
            else:
                rep1_paired_end = args.rep1pe
            if not simplicate_experiment:
                exp_rep2_ta = dxpy.dxlink(resolve_file(args.rep2[0]).get_id())
                exp_rep2_ta_desc = dxpy.describe(exp_rep2_ta)
                exp_rep2_mapping_analysis_id = dxpy.describe(exp_rep2_ta_desc['createdBy']['job'])['analysis']
                exp_rep2_mapping_analysis = dxpy.describe(exp_rep2_mapping_analysis_id)
                rep2_xcor_stage_description = next(
                    stage
                    for stage in exp_rep2_mapping_analysis.get('stages')
                    if stage['execution']['executableName'] == 'xcor')
                exp_rep2_cc = rep2_xcor_stage_description['execution']['output']['CC_scores_file']
                if args.rep2pe is None:
                    print("Inferring rep2 PE-ness from analysis")
                    rep2_paired_end = rep1_xcor_stage_description['execution']['output']['paired_end']
                else:
                    rep2_paired_end = args.rep1pe
            else:
                exp_rep2_ta = None
                exp_rep2_cc = None
                rep2_paired_end = None

            ctl_rep1_ta = dxpy.dxlink(resolve_file(args.ctl1[0]).get_id())
            if not unary_control and not simplicate_experiment:
                ctl_rep2_ta = dxpy.dxlink(resolve_file(args.ctl2[0]).get_id())
            else:
                ctl_rep2_ta = None
        else:  # blank workflow
            ctl_rep1_ta = None
            ctl_rep2_ta = None

            # here we need to calculate the cc scores files, because we're only
            # being supplied tagAligns
            # if we had mapped everything above we'd already have a handle to
            # the cc file
            xcor_only_applet = find_applet_by_name(XCOR_ONLY_APPLET_NAME, applet_project.get_id())
            # xcor_output_folder = resolve_folder(output_project, output_folder + '/' + xcor_only_applet.name)
            xcor_output_folder = xcor_only_applet.name
            xcor_only_stages = []
            rep1_xcor_input = {'spp_version': args.spp_version}
            if args.rep1pe is not None:
                rep1_xcor_input.update({'paired_end': args.rep1pe})
            exp_rep1_cc_stage_id = workflow.add_stage(
                xcor_only_applet,
                name="Rep1 cross-correlation",
                folder=xcor_output_folder,
                stage_input=rep1_xcor_input
            )
            xcor_only_stages.append({'xcor_only_rep1_id': exp_rep1_cc_stage_id})
            exp_rep1_cc = dxpy.dxlink(
                {'stage': exp_rep1_cc_stage_id,
                 'outputField': 'CC_scores_file'})
            rep1_paired_end = dxpy.dxlink(
                {'stage': exp_rep1_cc_stage_id,
                 'outputField': 'paired_end'})
            exp_rep1_ta = dxpy.dxlink(
                {'stage': exp_rep1_cc_stage_id,
                 'inputField': 'input_tagAlign'})
            if not simplicate_experiment:
                rep2_xcor_input = {'spp_version': args.spp_version}
                if args.rep2pe is not None:
                    rep2_xcor_input.update({'paired_end': args.rep2pe})
                exp_rep2_cc_stage_id = workflow.add_stage(
                    xcor_only_applet,
                    name="Rep2 cross-correlation",
                    folder=xcor_output_folder,
                    stage_input=rep2_xcor_input
                )
                xcor_only_stages.append({'xcor_only_rep2_id': exp_rep2_cc_stage_id})
                exp_rep2_cc = dxpy.dxlink(
                    {'stage': exp_rep2_cc_stage_id,
                     'outputField': 'CC_scores_file'})
                rep2_paired_end = dxpy.dxlink(
                    {'stage': exp_rep2_cc_stage_id,
                     'outputField': 'paired_end'})
                exp_rep2_ta = dxpy.dxlink(
                    {'stage': exp_rep2_cc_stage_id,
                     'inputField': 'input_tagAlign'})

            else:
                exp_rep2_cc = None
                exp_rep2_ta = None
                rep2_paired_end = None

    if not args.maponly:
        encode_macs2_applet = find_applet_by_name(ENCODE_MACS2_APPLET_NAME, applet_project.get_id())
        encode_macs2_stages = []
        # peaks_output_folder = resolve_folder(output_project, output_folder + '/' + encode_macs2_applet.name)
        peaks_output_folder = encode_macs2_applet.name

        # for simplicate experiments and/or unary controls, some of the ta inputs
        # will have the value None
        macs2_stage_input_mapping = {
                'rep1_ta' : exp_rep1_ta,
                'rep2_ta' : exp_rep2_ta,
                'ctl1_ta': ctl_rep1_ta,
                'ctl2_ta' : ctl_rep2_ta,
                'rep1_xcor' : exp_rep1_cc,
                'rep2_xcor' : exp_rep2_cc,
                'rep1_paired_end': rep1_paired_end,
                'rep2_paired_end': rep2_paired_end,
                'narrowpeak_as': dxpy.dxlink(resolve_file(args.narrowpeak_as)),
                'gappedpeak_as': dxpy.dxlink(resolve_file(args.gappedpeak_as)),
                'broadpeak_as':  dxpy.dxlink(resolve_file(args.broadpeak_as)),
                'genomesize': genomesize,
                'chrom_sizes': chrom_sizes,
                'fragment_length': args.fragment_length
            }

        # have to prune out any arguments with value None because DX will error
        # with arguments with null values
        macs2_stage_input = dict([(k,v) for k,v in macs2_stage_input_mapping.iteritems() if v is not None])

        encode_macs2_stage_id = workflow.add_stage(
            encode_macs2_applet,
            name='ENCODE Peaks',
            folder=peaks_output_folder,
            stage_input=macs2_stage_input
            )
        encode_macs2_stages.append({'name': 'ENCODE Peaks', 'stage_id': encode_macs2_stage_id})

        if run_idr:
            encode_spp_applet = find_applet_by_name(ENCODE_SPP_APPLET_NAME, applet_project.get_id())
            encode_spp_stages = []
            # idr_peaks_output_folder = resolve_folder(output_project, output_folder + '/' + encode_spp_applet.name)
            idr_peaks_output_folder = encode_spp_applet.name
            PEAKS_STAGE_NAME = 'SPP Peaks'
            # for simplicate experiments and/or unary controls, some of the ta inputs
            # will have the value None
            peaks_stage_input_mapping = {
                        'rep1_ta' : exp_rep1_ta,
                        'rep2_ta' : exp_rep2_ta,
                        'ctl1_ta': ctl_rep1_ta,
                        'ctl2_ta' : ctl_rep2_ta,
                        'rep1_xcor' : exp_rep1_cc,
                        'rep2_xcor' : exp_rep2_cc,
                        'rep1_paired_end': rep1_paired_end,
                        'rep2_paired_end': rep2_paired_end,
                        'as_file': dxpy.dxlink(resolve_file(args.narrowpeak_as)),
                        'idr_peaks': True,
                        'spp_version': args.spp_version,
                        'spp_instance': args.spp_instance,
                        'fragment_length': args.fragment_length
                        }
            if chrom_sizes:
                peaks_stage_input_mapping.update({'chrom_sizes': chrom_sizes})
            else:
                peaks_stage_input_mapping.update({'chrom_sizes': dxpy.dxlink({'stage': encode_macs2_stage_id, 'inputField': 'chrom_sizes'})})
            # have to prune out any arguments with value None because DX will error
            # with arguments with null values
            peaks_stage_input = dict([(k,v) for k,v in peaks_stage_input_mapping.iteritems() if v is not None])

            encode_spp_stage_id = workflow.add_stage(
                encode_spp_applet,
                name=PEAKS_STAGE_NAME,
                folder=idr_peaks_output_folder,
                stage_input=peaks_stage_input
                )
            encode_spp_stages.append({'name': PEAKS_STAGE_NAME, 'stage_id': encode_spp_stage_id})

            # TODO here I think we should abstract out all the IDR to one step like the two peak-calling steps
            idr_applet = find_applet_by_name(IDR2_APPLET_NAME, applet_project.get_id())
            encode_idr_applet = find_applet_by_name(ENCODE_IDR_APPLET_NAME, applet_project.get_id())
            idr_stages = []
            # idr_output_folder = resolve_folder(output_project, output_folder + '/' + idr_applet.name)
            idr_output_folder = idr_applet.name
            if (args.rep1 and args.ctl1 and args.rep2) or blank_workflow or simplicate_experiment:
                idr_stage_id = workflow.add_stage(
                    idr_applet,
                    name='IDR Rep 1 Self-pseudoreplicates',
                    folder=idr_output_folder,
                    stage_input={
                        'rep1_peaks' : dxpy.dxlink(
                            {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                             'outputField': 'rep1pr1_peaks'}),
                        'rep2_peaks' : dxpy.dxlink(
                            {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                             'outputField': 'rep1pr2_peaks'}),
                        'pooled_peaks': dxpy.dxlink(
                            {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                             'outputField': 'rep1_peaks'})
                    }
                )
                idr_stages.append({'name': 'IDR Rep 1 Self-pseudoreplicates', 'stage_id': idr_stage_id})

                if not simplicate_experiment:
                    idr_stage_id = workflow.add_stage(
                        idr_applet,
                        name='IDR Rep 2 Self-pseudoreplicates',
                        folder=idr_output_folder,
                        stage_input={
                            'rep1_peaks' : dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'rep2pr1_peaks'}),
                            'rep2_peaks' : dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'rep2pr2_peaks'}),
                            'pooled_peaks': dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'rep2_peaks'})
                        }
                    )
                    idr_stages.append({'name': 'IDR Rep 2 Self-pseudoreplicates', 'stage_id': idr_stage_id})

                    idr_stage_id = workflow.add_stage(
                        idr_applet,
                        name='IDR True Replicates',
                        folder=idr_output_folder,
                        stage_input={
                            'rep1_peaks' : dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'rep1_peaks'}),
                            'rep2_peaks' : dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'rep2_peaks'}),
                            'pooled_peaks': dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'pooled_peaks'})
                        }
                    )
                    idr_stages.append({'name': 'IDR True Replicates', 'stage_id': idr_stage_id})

                    idr_stage_id = workflow.add_stage(
                        idr_applet,
                        name='IDR Pooled Pseudoreplicates',
                        folder=idr_output_folder,
                        stage_input={
                            'rep1_peaks' : dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'pooledpr1_peaks'}),
                            'rep2_peaks' : dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'pooledpr2_peaks'}),
                            'pooled_peaks': dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_spp_stages if ss['name'] == PEAKS_STAGE_NAME),
                                 'outputField': 'pooled_peaks'})
                        }
                    )
                    idr_stages.append({'name': 'IDR Pooled Pseudoreplicates', 'stage_id': idr_stage_id})

                final_idr_stage_input = {
                        'r1pr_peaks': dxpy.dxlink(
                            {'stage': next(ss.get('stage_id') for ss in idr_stages if ss['name'] == 'IDR Rep 1 Self-pseudoreplicates'),
                             'outputField': 'IDR_peaks'}),
                        'rep1_ta': exp_rep1_ta,
                        'rep1_xcor': exp_rep1_cc,
                        'paired_end': rep1_paired_end,  # applies to replicated experiments, too
                        'as_file': dxpy.dxlink(resolve_file(args.narrowpeak_as)),
                        'rep1_signal': dxpy.dxlink(
                            {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == 'ENCODE Peaks'),
                             'outputField': 'rep1_fc_signal'})
                    }
                if not simplicate_experiment:
                    final_idr_stage_input.update({
                            'reps_peaks' : dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in idr_stages if ss['name'] == 'IDR True Replicates'),
                                 'outputField': 'IDR_peaks'}),
                            'r2pr_peaks' : dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in idr_stages if ss['name'] == 'IDR Rep 2 Self-pseudoreplicates'),
                                 'outputField': 'IDR_peaks'}),
                            'pooledpr_peaks': dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in idr_stages if ss['name'] == 'IDR Pooled Pseudoreplicates'),
                                 'outputField': 'IDR_peaks'}),
                            'rep2_ta': exp_rep2_ta,
                            'rep2_xcor': exp_rep2_cc,
                            'rep2_signal': dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == 'ENCODE Peaks'),
                                 'outputField': 'rep2_fc_signal'}),
                            'pooled_signal': dxpy.dxlink(
                                {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == 'ENCODE Peaks'),
                                 'outputField': 'pooled_fc_signal'})
                        })

                if blacklist:
                    final_idr_stage_input.update({'blacklist': blacklist})
                if chrom_sizes:
                    final_idr_stage_input.update({'chrom_sizes': chrom_sizes})
                else:
                    final_idr_stage_input.update({'chrom_sizes': dxpy.dxlink({'stage': encode_spp_stage_id, 'inputField': 'chrom_sizes'})})
                final_idr_stage_id = workflow.add_stage(
                    encode_idr_applet,
                    name='Final IDR peak calls',
                    folder=idr_output_folder,
                    stage_input=final_idr_stage_input,

                )
                idr_stages.append({'name': 'Final IDR peak calls', 'stage_id': final_idr_stage_id})

        if target_type == 'histone':
            PEAKS_STAGE_NAME = "ENCODE Peaks"
            overlap_peaks_applet = find_applet_by_name(OVERLAP_PEAKS_APPLET_NAME, applet_project.get_id())
            overlap_peaks_stages = []
            for peaktype in ['narrowpeaks', 'gappedpeaks', 'broadpeaks']:

                if peaktype == 'narrowpeaks':
                    as_file = dxpy.dxlink(resolve_file(args.narrowpeak_as))
                    peak_type_extension = 'narrowPeak'

                elif peaktype == 'gappedpeaks':
                    as_file = dxpy.dxlink(resolve_file(args.gappedpeak_as))
                    peak_type_extension = 'gappedPeak'

                elif peaktype == 'broadpeaks':
                    as_file = dxpy.dxlink(resolve_file(args.broadpeak_as))
                    peak_type_extension = 'broadPeak'

                overlap_peaks_stage_input = {
                    'rep1_peaks': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'rep1_%s' % (peaktype)}),
                    'rep2_peaks': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'rep2_%s' % (peaktype)}),
                    'pooled_peaks': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'pooled_%s' % (peaktype)}),
                    'pooledpr1_peaks': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'pooledpr1_%s' % (peaktype)}),
                    'pooledpr2_peaks': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'pooledpr2_%s' % (peaktype)}),
                    'rep1_ta': exp_rep1_ta,
                    'rep1_xcor': exp_rep1_cc,
                    'rep2_ta': exp_rep2_ta,
                    'rep2_xcor': exp_rep2_cc,
                    'paired_end': rep1_paired_end,  # applies to replicated experiments, too
                    'as_file': as_file,
                    'peak_type': peak_type_extension,
                    'prefix': 'final',
                    'rep1_signal': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'rep1_fc_signal'}),
                    'rep2_signal': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'rep2_fc_signal'}),
                    'pooled_signal': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'pooled_fc_signal'})
                } if not simplicate_experiment else {
                    'rep1_peaks': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'rep1pr1_%s' % (peaktype)}),
                    'rep2_peaks': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'rep1pr2_%s' % (peaktype)}),
                    'pooled_peaks': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'rep1_%s' % (peaktype)}),
                    'rep1_ta': exp_rep1_ta,
                    'rep1_xcor': exp_rep1_cc,
                    'paired_end': rep1_paired_end,  # applies to replicated experiments, too
                    'as_file': as_file,
                    'peak_type': peak_type_extension,
                    'prefix': 'final',
                    'rep1_signal': dxpy.dxlink(
                        {'stage': next(ss.get('stage_id') for ss in encode_macs2_stages if ss['name'] == PEAKS_STAGE_NAME),
                         'outputField': 'rep1_fc_signal'})
                }
                if chrom_sizes:
                    overlap_peaks_stage_input.update({'chrom_sizes': chrom_sizes})
                else:
                    overlap_peaks_stage_input.update({'chrom_sizes': dxpy.dxlink({'stage': encode_macs2_stage_id, 'inputField': 'chrom_sizes'})})

                overlap_peaks_stage_id = workflow.add_stage(
                    overlap_peaks_applet,
                    name='Final %s' % (peaktype),
                    folder=peaks_output_folder,
                    stage_input=overlap_peaks_stage_input
                )
                overlap_peaks_stages.append({'name': 'Final %s' %(peaktype), 'stage_id': overlap_peaks_stage_id})

    if args.yes:
        if args.debug:
            analysis = workflow.run({}, folder=output_folder, priority='high', debug={'debugOn': ['AppInternalError', 'AppError']}, delay_workspace_destruction=True, allow_ssh=['*'])
        else:
            analysis = workflow.run({}, folder=output_folder, priority='normal')

        analysis.set_properties({
            "target_type": target_type,
            "unreplicated_experiment": str(simplicate_experiment),
            "unary_control": str(unary_control)
        })
        print("Running %s as %s" % (analysis.name, analysis.get_id()))

        if args.accession:
            accession_analysis_applet = find_applet_by_name(ACCESSION_ANALYSIS_APPLET_NAME, applet_project.get_id())
            accession_output_folder = '/' + accession_analysis_applet.name
            accession_job_input = {
                'analysis_ids': [analysis.get_id()],
                'wait_on_files': []
            }
            if args.fqcheck is not None:
                accession_job_input.update({'fqcheck' : args.fqcheck})
            if args.skip_control is not None:
                accession_job_input.update({'skip_control' : args.skip_control})
            if args.force_patch is not None:
                accession_job_input.update({'force_patch': args.force_patch})
            # assert accession_stage_input['wait_on_files'], "ERROR: workflow has no wait_on_files defined, so --accession is not supported."
            time.sleep(5)
            max_retries = 10
            retries = max_retries
            while retries:
                try:
                    accession_job = accession_analysis_applet.run(
                        accession_job_input,
                        name='Accession %s' % (analysis.name),
                        folder=accession_output_folder,
                        depends_on=analysis.describe()['dependsOn']
                    )
                except Exception as e:
                    logging.error("%s launching auto-accession ... %d retries left" % (e, retries))
                    time.sleep(5)
                    retries -= 1
                    continue
                else:
                    logging.info("Auto-accession will run as %s %s" % (accession_job.name, accession_job.get_id()))
                    break
            else:
                logging.error("Auto-accession failed with %s" % ())


if __name__ == '__main__':
    main()
