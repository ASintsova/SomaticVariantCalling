import argparse
import subprocess
import shlex
import shutil
#from pathlib import Path
import yaml
import os
import sys

# Options


def get_args():
    parser = argparse.ArgumentParser("Somatic Variant Calling Pipeline\n")
    parser.add_argument("-a", "--analysis",
                        help="Options: preprocess, call, annotate", required=True)
    parser.add_argument("-c", "--config", nargs='?', type=str,
                        help="Config file", required=False)
    parser.add_argument('-j', '--j', help="Parallel param in snakemake",
                        type=int, default=1, required=False)
    parser.add_argument('-np', '--np', help="Dry run", action='store_true',
                        required=False)
    parser.add_argument("-local", "--local", help='Run locally or on cluster',
                        action='store_true', required=False)
    parser.add_argument('-unlock', '--unlock', action ='store_true', required=False)
    return parser


def build_image(analysis, local=True):
    analysis_to_image = {'preprocess': ['somvar'],
                         'call': ['somvar'],
                         'annotate': ['somvar'],
                         'build': ['somvar']}
    if local:
        for image in analysis_to_image[analysis]:
            if not os.path.isfile("singularity_images/{}.img".format(image)):
                print("No image found, building {} image".format(image))
                cmd = shlex.split(
                    "sudo singularity build singularity_images/{}.img singularity_images/{}_singularity".format(
                        image, image))
                subprocess.call(cmd)
            else:
                print("{} image found".format(image))
        return analysis_to_image[analysis]
    else:
        if analysis == 'build':
            print("Can't build on LeoMed")
            sys.exit()
        else:
            if not os.path.isfile("singularity_images/{}.img".format(analysis_to_image[analysis])):
                print("No image found, Exiting")
                sys.exit()
            else:
                return analysis_to_image[analysis]


def merge_config_and_settings(config):

    user_configs = yaml.load(open(config))#, Loader=yaml.FullLoader)
    subprocess.call(["mkdir", "-p", user_configs['outputDir']])
    new_config_path = user_configs['outputDir'] + "/" + user_configs['projectName'] + ".yaml"
    if os.path.isfile(new_config_path):
        return new_config_path
    else:
        with open(new_config_path, "w") as new_config:
            new_config.write(open("configs/settings.yaml", "r").read())
            new_config.write(open(config, "r").read())
        return new_config_path


def clean(config):
    user_configs = yaml.load(open(config))#, Loader=yaml.FullLoader)
    print("Removing {} ".format(user_configs['outputDir']))
    shutil.rmtree(user_configs['outputDir'])



def main(args):

    print(args)
    analysis = args.analysis

    image = "singularity_images/{}.img".format(build_image(analysis)[0], args.local)

    if analysis == 'build':
        return "Done"

    user_config = args.config
    j = args.j
    np = '-np' if args.np else ''
    unlock = '--unlock' if args.unlock else ''
    local = args.local
    new_config = merge_config_and_settings(user_config)

    if local:
        cmd_str = "singularity exec {} snakemake --configfile {} " \
                  "-j {} " \
                  "{} {} {}".format(image, new_config, j, np, analysis, unlock)
        print(cmd_str)
        subprocess.call(shlex.split(cmd_str))
    else:
        leomed = yaml.load(open(new_config))['leomed']#, Loader=yaml.FullLoader)['leomed']
        nodes = leomed['nodes']
        mem = leomed['mem']
        wallTime = leomed['walltime']
        dataDir = leomed['data']
        cmd_str1 = "bsub -n {} -R 'rusage[mem={}]' -W {} ".format(nodes, mem, wallTime)
        cmd_str2 = "singularity exec -H $HOME -B {}:/opt/data {} ".format(dataDir, image)
        # todo make sure local path to image works on leomed
        cmd_str3 = "snakemake --configfile {} -j {} {} {} {}".format(new_config, j, np, analysis, unlock)
        cmd_str = cmd_str1 + cmd_str2 + cmd_str3
        print(cmd_str)
        # todo test this on leomed
        subprocess.call(shlex.split(cmd_str))

if __name__ == "__main__":
    args = get_args().parse_args()
    if args.analysis == 'clean':
        clean(args.config)
    else:
        main(args)





