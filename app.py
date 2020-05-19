# Margo Raijmakers
# 18-05-20

from flask import Flask, render_template, request
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
import re

app = Flask(__name__)


@app.route('/')
def get_info():
    """Deze functie haalt de informatie op voor de webapplicatie

    :return: het template voor de webapplicatie
    """
    seq = request.args.get("seq", "")
    type = seq_type(seq)
    if type == "Deze sequentie is een DNA-sequentie":
        mrna_seq = transcription(seq)
        aa_seq = translation(seq)
        return render_template("Afvinkopdracht 4 HTML DNA.html", type=type,
                               mrna_seq=mrna_seq, aa_seq=aa_seq)
    elif type == "Deze sequentie is een aminozuursequentie":
        blast(seq)
        result = parse()
        return render_template("Afvinkopdracht 4 HTML aminozuur.html",
                               type=type, result=result)
    else:
        return render_template("Afvinkopdracht 4 HTML base.html",
                               type=type)


def seq_type(seq):
    """Deze functie checkt wat voor sequentie de opgegeven sequentie is

    :param seq: str - de opgegeven sequentie
    :return: str - "Deze sequentie is" + het type sequentie
    """
    if seq == "":
        return ""
    # DNA check
    if re.search("[^ATCGN]", seq.upper()):
        # RNA check
        if re.search("[^AUCGN]", seq.upper()):
            # aminozuur check
            if re.search("[^ARNDCFQEGHILKMPSTWYV]", seq.upper()):
                return \
                    "Deze sequentie is geen DNA-, RNA- of aminozuursequentie"
            else:
                return "Deze sequentie is een aminozuursequentie"
        else:
            return "Deze sequentie is een RNA-sequentie"
    else:
        return "Deze sequentie is een DNA-sequentie"


def transcription(seq):
    """Deze functie transcribeert de opgegeven mRNA-sequentie

    :param seq: str - de opgegeven sequentie
    :return: "mRNA-sequentie: " + de mRNA-sequentie
    """
    coding_dna = Seq(seq, IUPAC.unambiguous_dna)
    mrna = coding_dna.transcribe()
    return "mRNA-sequentie: " + mrna.upper()


def translation(seq):
    """Deze functie transleert de opgegeven mRNA-sequentie

    :param seq: str - de opgegeven sequentie
    :return: "Aminozuur-sequentie: " + de aminozuursequentie
    """
    coding_dna = Seq(seq, IUPAC.unambiguous_dna)
    aa = coding_dna.translate(table=1)
    return "Aminozuursequentie: " + aa


def blast(seq):
    """Deze functie blast de opgegeven aminozuursequentie

    :param seq: str - de opgegeven sequentie
    :return: een bestand met de blastresultaten
    """
    result_handle = NCBIWWW.qblast('blastp', 'nr', seq)
    with open("blast.xml", "w") as out_handle:
        out_handle.write(result_handle.read())


def parse():
    """"Deze functie haalt de nuttige informatie uit het blastbestand

    :return: het meest waarschijnlijke gen waarvan de aminozuursequentie
    afkomstig is
    """
    with open("my_blast.xml", "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)
        result = re.search(r'(?<=\| ).+?(?=\[)', str(blast_record.alignments[0]
                                                     .title)).group()
        return "Meest waarschijnlijke gen waar het eiwit van afkomstig is: " +\
               result


if __name__ == '__main__':
    app.run()
