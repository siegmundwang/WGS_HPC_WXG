aws s3 cp *_val_* *bam s3://jiguang2021/RNA_RESULT \
--endpoint-url=http://tos-s3-cn-qingdao.volces.com

aws s3 ls s3://jiguang2021/WGS --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com
aws s3 cp . s3://jiguang2021/RNA_RESULT/12TS0H0076NR3/ --recursive --exclude "12TS0H0076NR3_[12].fq.gz" --include "*.bam" \
    --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com


aws s3 rm s3://jiguang2021/RNA_RESULT/12TS0H0076NR3/200908_SEQ115_DP8400012879BR_L01_SP2008190773 --recursive \
    --endpoint-url=http://tos-s3-cn-qingdao.volces.com


rsem-calculate-expression --paired-end --no-bam-output --alignments -p 88 -q \
                12TS0H0076NR3_Aligned.toTranscriptome.out.bam /data/common_data/reference/RSEM_INDEX_GRCh38/GRCh38 RSEM/12TS0H0076NR3  


                14TS0H0110TR1

aws s3 cp  --quiet --recursive s3://jiguang2021/WGS/14TS0H0110TR1 /local_data/WGS/14TS0H0110TR1 \
        --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com # ivolces for 内网；volces for 外网
aws s3 cp  --quiet --recursive s3://jiguang2021/ReSequence/14TS0H0110TR1 /local_data/WGS/14TS0H0110TR1 \
        --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com # ivolces for 内网；volces for 外网

aws s3 cp nohup.out s3://jiguang2021/ --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com

# recursive 很重要！
aws s3 cp . s3://jiguang2021/DNA_RESULT/14TS0H0110TR1/ --recursive --exclude "*" --include "out_*_[12].fq.gz" --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com 


aws s3 cp . s3://jiguang2021/DNA_RESULT/ --exclude "*" --include "WGS*" --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com

aws s3api head-object --bucket jiguang2021 --key TEST_FILE --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com


aws s3api head-object --bucket jiguang2021 --key DNA_RESULT/14TS0H0110TR1/ --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com

yagmail -user xwang969-c@my.cityu.edu.hk -password qazwsxSd112112 -subject RNA-Seq_Condition \
    -contents `wc -l /data/Manager/RNA-Seq_Manage/Finished_STAR` \
    -to xwang969-c@my.cityu.edu.hk

usage: yagmail [-h] [-to TO [TO ...]] [-subject SUBJECT [SUBJECT ...]] [-contents CONTENTS [CONTENTS ...]]
               [-attachments ATTACHMENTS [ATTACHMENTS ...]] [-user USER] [-password PASSWORD]

Send a (g)mail with yagmail.

optional arguments:
  -h, --help            show this help message and exit
  -to TO [TO ...], -t TO [TO ...]
                        Send an email to address "TO"
  -subject SUBJECT [SUBJECT ...], -s SUBJECT [SUBJECT ...]
                        Subject of email
  -contents CONTENTS [CONTENTS ...], -c CONTENTS [CONTENTS ...]
                        Contents to send
  -attachments ATTACHMENTS [ATTACHMENTS ...], -a ATTACHMENTS [ATTACHMENTS ...]
                        Attachments to attach
  -user USER, -u USER   Username
  -password PASSWORD, -p PASSWORD
                        Preferable to use keyring rather than password here


    aws s3 cp . s3://jiguang2021/DNA_RESULT/${i}/ --recursive  --exclude "*"  --include "score.txt.idx" \
            --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com  && (echo "SUCCESS upload BAMs")

aws s3 cp  --quiet --recursive s3://jiguang2021/check . \
        --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com # ivolces for 内网；volces for 外网



        aws s3api head-object --bucket jiguang2021 --key WGS/14TSLC0127NR3/200904_SEQ103_DP8400012853TL_L01_SP2008070630/DP8400012853TL_L01_575_1.fq.gz \
             --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com

aws s3 cp  --quiet  s3://jiguang2021/WGS/${sample}/${sample}_Aligned.out.bam . \
        --endpoint-url=http://tos-s3-cn-qingdao.volces.com # ivolces for 内网；volces for 外网

16TSRC0664TD1/200810_SEQ124_DP8400012364TL_L01_SP1905200435/DP8400012364TL_L01_538_2.fq.gz


        aws s3api head-object --bucket jiguang2021 --key WGS/14TSRC0454NR3/200825_SEQ055_DP8400012596TR_L01_SP1910110512/DP8400012596TR_L01_521_1.fq.gz \
             --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com

aws s3 cp  --recursive s3://jiguang2021/ReSequence . --exclude "*" --include "*txt" \
        --endpoint-url=http://tos-s3-cn-qingdao.ivolces.com # ivolces for 内网；volces for 外网


STAR --runThreadN 128 --runMode genomeGenerate \
        --genomeFastaFiles .GRCh38.primary_assembly.genome.fa \
        --sjdbGTFfile Homo_sapiens.GRCh38.101.gtf \
        --sjdbOverhang 99 --genomeDir ./STAR_INDEX_GRCh38_2.7.8a/



featureCounts HKCI-C2_Aligned.sortedByCoord.out.bam -p -t exon -T 64 -o HKCI-C2_FC \
    -a /dataserver145/genomics/wangxg/common_data/reference/Homo_sapiens.GRCh38.101.gtf \
    -g gene_biotype



featureCounts \
    -B -C -g gene_biotype -t exon \
    -p \
    -T 64 \
    -a /dataserver145/genomics/wangxg/common_data/reference/Homo_sapiens.GRCh38.101.gtf \
    -s 0 \
    -o HKCI-C2_FC \
    HKCI-C2_Aligned.sortedByCoord.out.bam

cut -f 1,7 HKCI-C2_FC | tail -n +3 | cat biotypes_header.txt - >> HKCI-C2.biotype_counts_mqc.tsv
mqc_features_stat.py HKCI-C2.biotype_counts_mqc.tsv -s HKCI-C2 -f rRNA -o HKCI-C2.biotype_counts_rrna_mqc.tsv



# install_kubes.sh
# 无包管理
CNI_VERSION="v0.8.2"
apt install ebtables ethtool socat conntrack
sudo mkdir -p /opt/cni/bin
curl -L "https://github.com/containernetworking/plugins/releases/download/${CNI_VERSION}/cni-plugins-linux-amd64-${CNI_VERSION}.tgz" | sudo tar -C /opt/cni/bin -xz
DOWNLOAD_DIR=/usr/local/bin
sudo mkdir -p $DOWNLOAD_DIR
CRICTL_VERSION="v1.17.0"
curl -L "https://github.com/kubernetes-sigs/cri-tools/releases/download/${CRICTL_VERSION}/crictl-${CRICTL_VERSION}-linux-amd64.tar.gz" | sudo tar -C $DOWNLOAD_DIR -xz

RELEASE="$(curl -sSL https://dl.k8s.io/release/stable.txt)"
cd $DOWNLOAD_DIR
sudo curl -L --remote-name-all https://storage.googleapis.com/kubernetes-release/release/${RELEASE}/bin/linux/amd64/{kubeadm,kubelet,kubectl}
sudo chmod +x {kubeadm,kubelet,kubectl}

RELEASE_VERSION="v0.4.0"
curl -sSL "https://raw.githubusercontent.com/kubernetes/release/${RELEASE_VERSION}/cmd/kubepkg/templates/latest/deb/kubelet/lib/systemd/system/kubelet.service" | sed "s:/usr/bin:${DOWNLOAD_DIR}:g" | sudo tee /etc/systemd/system/kubelet.service
sudo mkdir -p /etc/systemd/system/kubelet.service.d
curl -sSL "https://raw.githubusercontent.com/kubernetes/release/${RELEASE_VERSION}/cmd/kubepkg/templates/latest/deb/kubeadm/10-kubeadm.conf" | sed "s:/usr/bin:${DOWNLOAD_DIR}:g" | sudo tee /etc/systemd/system/kubelet.service.d/10-kubeadm.conf


systemctl enable --now kubelet


apt install ebtables ethtool socat conntrack


aws s3 cp   s3://jiguang2021/RNA_RESULT/12TS0H0076NR3/12TS0H0076NR3_Aligned.out.bam . \
  --endpoint-url=http://tos-s3-cn-qingdao.volces.com # ivolces for 内网；volces for 外网


  aws s3 cp  s3://jiguang2021/RNA_RESULT/12TS0H0076NR3/12TS0H0076NR3_Aligned.out.bam 12TS0H0076NR3_Aligned.out.bam --endpoint-url=http://tos-s3-cn-qingdao.volces.com



 featureCounts -B -C -g gene_biotype -t exon -p -T 64 \
        -a /data/common_data/reference/Homo_sapiens.GRCh38.101.gtf \
        -s 0 -o 12TS0H0076TR3_FC 12TS0H0076TR3_Aligned.out.bam
