import subprocess
import re
import requests
import yaml

with open("config.yaml", "r") as yamlfile:
    data = yaml.load(yamlfile, Loader=yaml.FullLoader)
    myvirus = data['virusname']


def downloader(family):
    """ Downloads the HTML page for the virus family of choice"""
    # Gets the correct URL and Virus family from the config
    url = f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name={family}"
    geturl = requests.get(url)

    # Writes the entire HTML page into a txt format
    with open(f"{family}.txt", "w") as out_f:
        out_f.write(geturl.text.strip())


def getter(family):
    """ Reworks the HTML txt file into all the viruses belonging to the chosen family"""
    viruses = []
    with open(f"{family}.txt") as htmlfile:
        for line in htmlfile:
            # All viruses in NCBI have links so here we filter all the links.
            if line.startswith("<LI"):
                # Almost all links have the same start and before line 115 ~ 120 no names are mentioned
                first = line.strip()[115:].split(">", )[1]
                # Uses a regular expression to remove anything behind the virus name such as (Arizona/2019)
                second = re.sub("[\(\[].*?[\)\]]+", "", first.split("<")[0])
                # Append all the viruses to a list for later writing.
                viruses.append(" ".join(second.split()[:8]))
                # Clears the list from all duplicates.
                viruses = list(dict.fromkeys(viruses))

    # Opens a new file for the writing of all the virusnames
    with open("Keywords.txt", "w") as out_file:
        for index, line in enumerate(viruses):
            # Makes sure to skip empty lines that .strip does not find.
            if len(line) >= 2:
                # Adds an enter after every virus so they dont clutter together.
                out_file.write(line.strip() + "\n")
        # Closes the file as per standard regulations.
        out_file.close()
        htmlfile.close()
    # Removes the HTML txt file since it is now obsolete
    subprocess.call(f"rm {family}.txt", shell=True)


def main():
    downloader(myvirus)
    getter(myvirus)
    yamlfile.close()


# Calls upon the main to make the functions work.
if __name__ == "__main__":
    main()
