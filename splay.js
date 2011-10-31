var fs    = require("fs");
var cl    = require("termcolor").define();
var LS    = require("linestream");
var AP    = require('argparser');
var Junjo = require("junjo");
var spawn = require("child_process").spawn;

// EXON FLAGS
const SELF_FLAG     = 0x01;  // in the same splicing ID
const GROUP_FLAG    = 0x02;  // in the same splicing GROUP
const KNOWNSET_FLAG = 0x04;  // known overlap
const TERM_FLAG     = 0x08;  // has terminate signal
const SELFTERM_FLAG = 0x10;  // same splicing ID with terminate signal

// function for flag checking
function has(flag, val) {
  return !! (flag & val);
}

/**
 * EXON BED INFORMATION
 *
 * <BASIC INFO>
 * 0  rname
 * 1  start
 * 2  end
 * 3  strand
 * 4  ID
 * 5  the number of patterns being contained
 * 6  the number of the first exon patterns
 * 7  the number of the last exon patterns
 * 8  is protein or not
 * 9  exon position information of each splicing pattern
 * 10 protein ID or *
 * --------------------------------------
 * <SV INFO>
 * 11 fragment id
 * 12 fragment total 
 * 13 SV type
 * 14 original rname
 * 15 original start
 * 16 original end
 * 17 original strand
 *
 * <EXAMPLE>
 * chr1	19941	20036	-	33 2	0	0	0	uc009viy.2[4/9],uc009vjf.2[5/7]	*	1	1	DUP	chr1	18267	18362	-
 **/

function main() {
  const p = new AP().addOptions([]).addValueOptions([]).parse();
  var convertedFile = p.getArgs(0);
  var originalFile  = p.getArgs(1);

  function showUsage() {
    const cmd = p.getOptions('exename') || (process.argv[0] + ' ' + require('path').basename(process.argv[1]));
    console.error('[synopsis]');
    console.egreen('\t' + cmd + ' <converted exon bed file> <original exon bed file>');
    console.error('[options]');
    console.error('');
  }

  if (!require("path").existsSync(convertedFile)) {
    if (!convertedFile) console.error("please specify a converted exon bed file.");
    else console.error(convertedFile, "no such file");
    return showUsage();
  }

  if (!require("path").existsSync(originalFile)) {
    if (!originalFile) console.error("please specify a original exon bed file.");
    else console.error(originalFile, "no such file");
    return showUsage();
  }

  var $j = new Junjo({noTimeout: true});

  $j.inputs(["convertedFile", "originalFile"]);

  $j("exonInfo", function(filename) {
    getExonInfo(filename, this.cb);
  })
  .using("originalFile");

  $j("exons", function(filename) {
    return fs.readFile(filename, "utf-8", this.cb);
  })
  .post(function(data) {
    return data.split("\n").map(function(v) {
      return v.split("\t");
    })
    .filter(function(data) {
      if (data.length < 11) return false;
      if (data.length < 12) return true;

      // filter incomplete fragments by SV type
      return (["ADUP", "ATRA", "*"].indexOf(data[13]) >= 0);
    });
  })
  .using("convertedFile")
  .eshift();

  $j("check", function(exons, groups, rels) {
    this.$.seqs = {};

    // for every exon
    this.forEach(exons, function(data, k) {
      var rname  = data[0];
      var strand = data[3];
      var exonID = data[4];
      var delta = (strand == "+") ? 1 : -1;
      var startCount = Number(data[6]);
      if (!startCount) return;

      // get start exons
      var starts = data[9].split(",").filter(function(v) {
        return v.match(/\[1\//);
      })
      .map(function(v) {
        return v.split("[")[0];
      });

      var spos = data[1];

      // for each start splice pattern
      starts.forEach(function(name, k2) {
        var id = name + "_" + exonID;
        var resultExons = {};
        var current = k - delta;

        // get exons
        while (true) {
          current += delta;

          // exons must exist
          var ex = exons[current];
          if (!ex) {
            var direc = (delta == 1) ? "below" : "above";
            console.ered("no more exons", direc, name, exonID, current);
            break;
          }

          // must be the same chromosome
          var thisRname = ex[0];
          if (rname!= thisRname) continue;

           // must be the same strand
          var thisStrand = ex[3];
          if (strand != thisStrand) continue;

          // exons must be within 2800000 base distance
          if (Math.abs(spos - ex[2]) > 2800000) {
            console.ered("exons must be within 2800000 base distance", Math.abs(spos - ex[2]), name, exonID, current);
            break;
          }

          // splice patterns in this exon
          var splices = ex[9].split(",");

          // put empty flag to current exon
          resultExons[current] = 0;

          // for each splice pattern in this exon, check termination flag
          splices.forEach(function(v) {
            var result = v.match(/\[([0-9]+)\/([0-9]+)\]/);
            var thisname = v.split("[")[0];
            if (!result || result[1] != result[2]) return;
            resultExons[current] |= TERM_FLAG; // set termination flag

            if (thisname == name)
              resultExons[current] |= SELFTERM_FLAG; // set selfterm flag
          });

          // list of splice pattern names
          var splicenames = splices.map(function(v) {
            return v.split("[")[0];
          });
          if (!splicenames.length) continue;

          // if containing self, set flag
          if (splicenames.indexOf(name) >= 0) {
            resultExons[current] |= SELF_FLAG;

            // if this exon has the terminate signal of the target splicing pattern, finish fetching
            if (has(resultExons[current], SELFTERM_FLAG)) {
              break;
            }
            continue; // no need to check same group or known relation
          }

          // checking group, relation
          splicenames.forEach(function(thisname, i) {
            var key = [name, thisname].sort().join(",");
            if (groups[key]) resultExons[current] |= GROUP_FLAG;
            if (rels[key]) resultExons[current] |= KNOWNSET_FLAG;
          }, this);
        }
        this.$.seqs[id] = resultExons;
      }, this);
    });
  })
  .after("exons", "exonInfo");

  $j.on("end", function(err, out) {
    var exons = out.exons;

    // for each sequences
    Object.keys(this.$.seqs).forEach(function(id) {
      var exs = this.$.seqs[id];
      var k = 0;

      // for each exons in the sequence until termination
      Object.keys(exs).every(function(exon_num) {
        var flag = exs[exon_num];

        // known relation
        if (has(flag, KNOWNSET_FLAG)) {
          // console.eblue("known relation, skipped", exon_num);
          return true;
        }

        // same group nonself exon
        if (has(flag, GROUP_FLAG) && !has(flag, SELF_FLAG)) {
          // console.eblue("same group, skipped", exon_num);
          return true;
        }

        // here remains UNKNOWN RELS or SELF. display it.
        var exon = exons[exon_num];
        console.log(exon.slice(0,5).concat([id, ++k, flag]).join("\t"));

        // if termination
        if (has(flag, SELFTERM_FLAG) || ( has(flag, TERM_FLAG) && !has(flag, SELF_FLAG))) {
          return false;
        }
        return true;
      }, this);
    }, this);

    console.egreen("finished");
  });

  $j.run(convertedFile, originalFile);
}


// get splicing groups, related splicing patterns
function getExonInfo(filename, callback) {
  var groups = {};
  var rels = {};

  var ex2gr = spawn("node", ["exon2groups.js", filename]);
  var ex2rl = spawn("node", ["exon2rels.js", filename]);

  ex2gr.stdout.setEncoding("utf8");
  ex2rl.stdout.setEncoding("utf8");

  var groupLines = new LS(ex2gr.stdout, {trim : true});

  groupLines.on("data", function(line) {
    var data = line.split("\t");
    if (data.length < 2) return;
    
    for (var i=0, l=data.length; i<l-1;i++) {
      for (var j=0; j<l; j++) {
        var key = [data[i], data[j]].sort().join(",");
        groups[key] = true;
        ex2rl.stdin.write(key + "\n");
      }
    }
  });

  groupLines.on("end", function() {
    ex2rl.stdin.end();
  });



  ex2rl.stderr.on("data", function(da) {
    console.log(da.toString());
  });

  var relLines = new LS(ex2rl.stdout, {trim : true});

  relLines.on("data", function(line) {
    var key = line.split("\t").sort().join(",");
    rels[key] = true;
  });

  relLines.on("end", function() {
    console.egreen("groups", Object.keys(groups).length);
    console.egreen("rels", Object.keys(rels).length);
    callback(groups, rels);
  });
}




if (process.argv[1] === __filename) { main(); }
