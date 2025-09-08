import argparse
import csv
import os
import pandas as pd
import urllib.request
import gzip
import time
import re
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--url_csv', type=str, help='file of GSM/GSE and URL')
    parser.add_argument('-f', '--folder_path', type=str, default = 'load_data',help='path to store data')
    parser.add_argument('-n', '--name_of_URL_column', type=str,default="GSE_URL",required=False, help='The column name of the column in which the URL is located in the csv file')
    parser.add_argument('-ne', '--name_of_GSE_column', type=str,default="GSE",required=False, help='The column name of the column in which the GSE is located in the csv file')
    args = parser.parse_args()
    url_csv = args.url_csv
    folder_path = args.folder_path
    name_of_URL_column = args.name_of_URL_column
    name_of_GSE_column = args.name_of_GSE_column
    table = pd.read_csv(url_csv, sep=",")
    series_list = table[name_of_GSE_column].tolist()
    result = table[table[name_of_GSE_column].isin(series_list)][[name_of_GSE_column, name_of_URL_column]]

    for i in result.index:
        start = time.time()
        gseid = result.loc[i, name_of_GSE_column]
        gse_three = gseid[0:6]
        path = "%s/%s/%s" % (folder_path, gse_three, gseid)
        print(path)
        url_list = result.loc[i, name_of_URL_column]
        url_list = url_list.strip("[]")
        url_list = url_list.split(", ")
        url_group = []
        for url_str in url_list:
            if isinstance(url_str, str):
                cleaned_url = url_str.strip("'")
                url_group.append(cleaned_url)  # 将每个字符串添加到 url_group 列表中
            else:
                continue
        print(url_group)
        try:
            os.makedirs(path)
        except:
            print("%s has created!" % path)
        for url in url_group:
            filename = url.split("/")[-1]
            if not os.path.exists("%s/%s" % (path, filename)):
                os.system("wget %s -O %s/%s" % (url, path, filename))
                print("finish download")
            else:
                print("file existed!")
        # download series matrix
        id_cut = gseid.__len__() - 3
        series_matrix_link = "https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/matrix/%s_series_matrix.txt.gz" % (
            gseid[0:id_cut], gseid, gseid)
        series_matrix_name = series_matrix_link.split("/")[-1]
        if not os.path.exists("%s/%s" % (path, series_matrix_name)):
            os.system("wget %s -O %s/%s" %
                    (series_matrix_link, path, series_matrix_name))
        try:
            filesize = os.path.getsize("%s/%s" % (path, series_matrix_name))
            if filesize == 0:
                os.system("rm %s/%s" % (path, series_matrix_link.split("/")[-1]))
                page_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/matrix/' % (
                    gseid[0:id_cut], gseid)
                html = urllib.request.urlopen(page_url).read().decode('utf-8')
                link_group = re.compile(
                    r'<a href=".*series_matrix\.txt\.gz"').findall(html)
                for link in link_group:
                    series_matrix_name = link.rstrip('"').lstrip('<a href="')
                    os.system("wget %s%s -O %s/%s" %
                            (page_url, series_matrix_name, path, series_matrix_name))
        except:
            # in case no file
            time.sleep(10)
            page_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/matrix/' % (
                gseid[0:id_cut], gseid)
            html = urllib.request.urlopen(page_url).read().decode('utf-8')
            link_group = re.compile(
                r'<a href=".*series_matrix\.txt\.gz"').findall(html)
            for link in link_group:
                series_matrix_name = link.rstrip('"').lstrip('<a href="')
                os.system("wget %s%s -O %s/%s" %
                        (page_url, series_matrix_name, path, series_matrix_name))
        with open('%s/genome_build.txt' % folder_path, 'a+') as build_file:
            with gzip.open('%s/%s' % (path, series_matrix_name), 'rb') as f:
                file_content = f.read()
                if b'hg38' in file_content or b'GRCh38' in file_content:
                    build_version = 'hg38'
                if b'hg19' in file_content or b'GRCh37' in file_content:
                    build_version = 'hg19'
                else: build_version = ''
                build_file.write('%s\t%s\n' % (gseid, build_version))
        end = time.time()
        print('Total time:', end - start)

if __name__ == "__main__":
    main()
    