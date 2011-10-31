#!/usr/bin/env node
var dna         = require("dna");
var FASTAReader = require("fastareader");
var SortedList = require("sortedlist");
var fs = require("fs");
var LS = require("linestream");
var AP = require("argparser");
var cl = require("termcolor").define();

const FASTA_HG19 = __dirname + "/../../data/hg19.clean.fasta";
const FASTA_JSON = FASTA_HG19 + ".json";
const FASTA_LINELEN = 50;


function main() {
  var p = new AP().parse();
  var splice_info_file = p.getArgs(0);
  var fastaFile = p.getArgs(1);
  var jsonFile  = p.getArgs(2);

  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + require('path').basename(process.argv[1]));
    console.error('[synopsis]');
    console.egreen('\t' + cmd + ' <splicing bed file> <fasta file> [<json file>]');
    console.error('[options]');
    console.error('');
  }


  if (!splice_info_file) {
    console.ered("requires splice info bed file");
    return showUsage();
  }

  if (!require("path").existsSync(splice_info_file)) {
    console.ered(splice_info_file, "no such file.");
    return showUsage();
  }

  if (!fastaFile) {
    console.ered("requires fasta file");
    return showUsage();
  }

  if (!require("path").existsSync(fastaFile)) {
    console.ered(fastaFile, "no such file.");
    return showUsage();
  }

  const filter = function(val, pos) {
    return (this.arr[pos]   == null || (this.arr[pos] != null && this.arr[pos][1] < val[0])) 
    &&     (this.arr[pos+1] == null || (this.arr[pos+1] != null && val[1] < this.arr[pos+1][0]));
  };
  const compare = function(a, b) {
    if (a == null) return -1;
    if (b == null) return  1;
    var c = a[0] - b[0];
    return (c > 0) ? 1 : (c == 0)  ? 0 : -1;
  };


  var lines = new LS(splice_info_file, {trim: true});

  var json = (function() {
    try {
      return JSON.parse(fs.readFileSync(jsonFile));
    } catch(e) {
      return null;
    }
  })();

  var freader = new FASTAReader(fastaFile, json);

  var currentList = null;
  var currentData = {};

  // chr1	12612	12721	+	2_1	uc001aaa.3_1_1	2	0
  lines.on("data", function(line) {
    if (!line || line.charAt(0) == "#") return;
    var data = line.split("\t");
    if (data.length < 5) throw new Error(data);

    var rname  = data[0];
    var start  = data[1];
    var end    = data[2];
    var strand = data[3];
    //var exonID = data[4];
    var splID  = data[5];
    //var order  = data[6];
    //var flag   = data[7];
    //var length = end - start + 1;

    if (currentData.ID != splID) {
      if (currentList) showData(currentList, currentData, freader);

      currentList = new SortedList(null, {filter: filter, compare: compare});
      process.stdout.write(">" + splID + "\n");
      currentData = {
        rname : rname,
        strand : strand,
        ID     : splID
      };
    }
    currentList.insert([start, end]);
  });

  lines.on("end", function() {
    showData(currentList, currentData, freader);
    freader.close();
  });

  function showData(list, data, reader) {
    try {
      var seq = list.toArray().reduce(function(sofar, val) {
        return sofar + reader.fetch(currentData.rname, val[0], val[1] - val[0] + 1);
      }, "");

      if (data.strand == "-") {
        // console.eyellow(data.ID, "rev");
        seq = dna.complStrand(seq, true);
      }
      // console.epurple("seq length", seq.length);
      while (seq) {
        console.log(seq.slice(0,FASTA_LINELEN));
        seq = seq.slice(FASTA_LINELEN);
      }
      console.log("");
    }
    catch (e) {
      console.ered(e.message);
    }
  }
}


if (process.argv[1].match('/([^/]+?)(\.js)?$')[1] == __filename.match('/([^/]+?)(\.js)?$')[1]) main();
