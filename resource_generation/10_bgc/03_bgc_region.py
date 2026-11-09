from Bio import SeqIO
import pandas as pd
import os

for i in range(8,10):
    folder_path = f'~/soil/15_bgc/cpsnj0{i}/all'

    for filename in os.listdir(folder_path):
        if filename.endswith(".gbk"):
            file_path = os.path.join(folder_path, filename)
            genome = filename.split('.')[0]
            output_file = genome + '_out_region.tsv'
            out_path = os.path.join(folder_path, output_file)

            # Collect all protoclusters and all regions from antiSMASH annotations
            all_protoclusters=[]
            all_regions=[]
            contig_n =1

            genbank = SeqIO.parse(file_path,format = "genbank")

            # get information for regions and protoclusters annotated by antiSMASH
            for record in genbank:

                # get information from regions annotated by antiSMASH
                regions = [f for f in record.features if f.type == "region" and f.qualifiers["tool"]==['antismash']]
                for region in regions:
                    r_num = "Region number " + str(contig_n)+ "." +  str(region.qualifiers["region_number"])[2:-1]
                    r_product = str(region.qualifiers["product"])[2:-1].replace("'","")
                    r_st=int(region.location.start)+1 #because indexing from 0 in python
                    r_end=int(region.location.end)
                    r_contig = record.id
                    r_genome = genome
                    r_info=[r_num, r_product, r_st, r_end, r_contig, r_genome]
                    all_regions.append(r_info)

                    # get information from protoclusters annotated by antiSMASH
                protoclusters = [f for f in record.features if f.type == "protocluster" and f.qualifiers["tool"]==['antismash']]
                for cluster in protoclusters:
                    c_num = "Protocluster number " + str(contig_n)+ "." +  str(cluster.qualifiers["protocluster_number"])[2:-1]
                    r_product = str(cluster.qualifiers["product"])[2:-1].replace("'","")
                    r_catogory = str(cluster.qualifiers["category"])[2:-1].replace("'","")
                    r_st=int(cluster.location.start)+1 #because indexing from 0 in python
                    r_end=int(cluster.location.end)
                    r_contig = record.id
                    r_info=[c_num, r_product,r_st, r_end, r_contig,r_catogory]
                    all_protoclusters.append(r_info)
                contig_n+=1

            # match protoclusters to regions
            all_total=[]

            for region in all_regions:
                r_num = region[0]
                r_product = region[1]
                r_start = region[2]
                r_stop= region[3]
                r_contig = region[4]
                r_genome = region[5]

                for proto in all_protoclusters:
                    c_start = proto[2]
                    c_stop = proto[3]
                    c_contig = proto[4]
                    c_catogory = proto[5]

                    if c_start >= r_start and c_stop <= r_stop and r_contig==c_contig:
                        new_row = proto
                        new_row.append(r_num)
                        new_row.append(r_product)
                        new_row.append(r_start)
                        new_row.append(r_stop)
                        new_row.append(r_contig)
                        new_row.append(r_genome)

                        all_total.append(new_row)

            # construct dataframe
            total_df = pd.DataFrame(all_total,columns=["protocluster","product","gbk_start","gbk_end","contig name","catogory","region","product","gbk_start","gbk_end","contig name",'genome'])

            # write collected dataframe to file
            total_df.to_csv(out_path)
