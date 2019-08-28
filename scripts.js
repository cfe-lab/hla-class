var interval = 0;
var done = false;
var error = false;
var seqWidth = 60;
var seqBlockSize = 10;
var batchMode = false;
var update = false;

window.onunload = function () {
	var request = new XMLHttpRequest();
	request.open("GET", "cleanup.php", false);
	request.send(null);
	return true;
}

function removeChildren(node) {
	while (node.hasChildNodes()) {
		node.removeChild(node.firstChild);
	}
}

function resetFile(id) {
	document.getElementById(id).setAttribute("type", "input");
	document.getElementById(id).setAttribute("type", "file");
}

function resetForm() {
	document.getElementById("fasta_text").value = "";
	document.getElementById("fullseq").value = "";
	document.getElementById("exon2").value = "";
	document.getElementById("exon3").value = "";

	resetFile("fasta_file");
	resetFile("fullseq_file");
	resetFile("exon2_file");
	resetFile("exon3_file");

	var request = new XMLHttpRequest();
	request.open("GET", "cleanup.php", false);

	clearMessages();
}

function updateAlleles() {
	clearMessages();
	done = false;
	error = false;

	var request = new XMLHttpRequest();
	request.open("GET", "cleanup.php", false);

	request.onreadystatechange = function () {
		if (request.readyState == 4 && request.status == 200) {
			document.getElementById("status").innerHTML="processing...";
			var files = request.responseText.split("\n");
			var updatesFile = files[0];
			var techErrorsFile = files[1];
			done = false;
			error = false;
			interval = setInterval(function() {
				getOutput(techErrorsFile, "tech_error", "tech_error", parseErrorsLine,
				"tech_error_header"); 
				getOutput(updatesFile, "update", "update", parseUpdateLine,
				"update_header");
				if (done) {
					clearInterval(interval);
					if (!error) {
						document.getElementById("status").innerHTML="done";
					} else {
						document.getElementById("status").innerHTML="error";	
					}
				}
			}, 200);
		}
	}
	request.open("POST", "update.php", true);
	request.send(null);
}

function switchMode(mode) {
	clearMessages();
	var batchTab = document.getElementById("batch_tab");
	var singleTab = document.getElementById("single_tab");
	var updateTab = document.getElementById("update_tab");
	var batchModeSpan = document.getElementById("batch_mode");
	var singleModeSpan = document.getElementById("single_mode");
	var inputSpan = document.getElementById("input");
	var submit = document.getElementById("submit");
	var updateModeSpan = document.getElementById("update_mode");

	if (mode == "batch") {
		batchMode = true;
		update = false;
		batchModeSpan.className = "";
		singleModeSpan.className = "hidden";
		inputSpan.className = "";
		updateModeSpan.className = "hidden";
		batchTab.className = "tab_selected";
		singleTab.className = "";
		updateTab.className = "";
		submit.onclick = sendForm;
	} else if (mode == "single") {
		batchMode = false;
		update = false;
		batchModeSpan.className = "hidden";
		singleModeSpan.className = "";
		updateModeSpan.className = "hidden";
		inputSpan.className = "";
		batchTab.className = "";
		singleTab.className = "tab_selected";
		updateTab.className = "";
		submit.onclick = sendForm;
	} else if (mode == "update") {
		update = true;
		updateModeSpan.className = "";
		batchModeSpan.className = "hidden";
		singleModeSpan.className = "hidden";
		inputSpan.className = "hidden";
		batchTab.className = "";
		singleTab.className = "";
		updateTab.className = "tab_selected";
		submit.onclick = updateAlleles;
	}
}

function validateMismatches(value) {
	var errnode = document.getElementById("form_errors");
	if (isNaN(parseInt(value)) || parseInt(value) < 0) {
		errnode.innerHTML = "Threshold must be a positive integer.";
	} else {
		while (errnode.hasChildNodes()) {
			errnode.removeChild(errnode.firstChild);
		}
		errnode.innerHTML = "";
	}
}

function getInput () {
	if (batchMode) {
		return document.getElementById("fasta_text").value.replace(/ /g, '');
	} else {
		var input_type = document.getElementById("input_type").value;
		if (input_type == "single") {
			var short = document.getElementById("fullseq_short").checked;
			var fullseq = document.getElementById("fullseq").value.replace(/ /g, '');
			if (fullseq.indexOf(">") == 0) {
				return fullseq;
			} 
			if (short) {
				return ">_short\n" + fullseq;
			} else {
				return ">_\n" + fullseq;
			}
		} else {
			var exon2_short = document.getElementById("exon2_short").checked;
			var input = "";
			var exon2 = document.getElementById("exon2").value.replace(/ /g, '');

			if (exon2.indexOf(">") != 0) {
				if (exon2_short) {
					input += ">_exon2_short\n";
				} else {
					input += ">_exon2\n";
				}
			}
			input += exon2;

			var exon3_short = document.getElementById("exon3_short").checked;
			var exon3 = document.getElementById("exon3").value.replace(/ /g, '');

			if (exon3.indexOf(">") != 0) {
				if (exon3_short) {
					input += "\n>_exon3_short";
				} else {
					input += "\n>_exon3";
				}
			}
			input += "\n";
			input += exon3;
			return input;
		}
	}
}

function toggleOptions (select) {
	var letter = select.value;
	var inputType = document.getElementById("input_type");
	if (letter == "A") {
		document.getElementById("input_type_choice").className = "hidden";
		inputType.value = "single";
		toggleInputType();
	} else {
		document.getElementById("input_type_choice").className = "";
	}
}

function toggleInputType(type) {
	var single = document.getElementById("single_mode_single");
	var linked = document.getElementById("single_mode_linked");

	if (type == "linked") {
		linked.className = "";
		single.className = "hidden";
	} else {
		single.className = "";
		linked.className = "hidden";
	}
}

function validateForm() {
	var errnode = document.getElementById("form_errors");
	while (errnode.hasChildNodes()) {
		errnode.removeChild(errnode.firstChild);
	}

	var nfiles = document.getElementById("fasta_file").files.length;
	var fasta_text = getInput();
	var ptn=/((^|\n)>.*\n([atcgyrwskmdvhbn-]+)(\n[atcgyrwskmdvhbn-]+)*)+/i;
	if (!ptn.test(fasta_text)) {
		var error = "Invalid sequences.";
		var node = document.createTextNode(error);
		errnode.appendChild(node);
		return false;
	}

	return true;
}

function clearMessages() {
	document.getElementById("status").innerHTML="";
	var divids = new Array("result", "error", "tech_error", "detail");
	var headers = new Array("result_header", "error_header", "tech_error_header",
		"detail_header");
	for (var i=0; i<4; ++i) {
		var div = document.getElementById(divids[i]);
		while (div.hasChildNodes()) {
			div.removeChild(div.firstChild);
		}
		var h = document.getElementById(headers[i]);
		h.className = "hidden";
	}
	document.getElementById("download").className = "hidden";
	var errnode = document.getElementById("form_errors");
	while (errnode.hasChildNodes()) {
		errnode.removeChild(errnode.firstChild);
	}
}

function sendForm() {
	if (!validateForm()) { return false; }
	clearMessages();
	done = false;
	error = false;

	var request = new XMLHttpRequest();
	var formData = new FormData();

	formData.append("fasta_text", getInput());
	formData.append("letter", document.getElementById("letter").value);

	if (!batchMode) {
		formData.append("threshold", document.getElementById("threshold").value);
	} else {
		formData.append("threshold", -1);
	}

	request.open("GET", "cleanup.php", false);

	request.onreadystatechange = function () {
		if (request.readyState == 4 && request.status == 200) {
			document.getElementById("status").innerHTML="processing...";
			var files = request.responseText.split("\n");
			var resultsFile = files[0];
			var detailsFile = files[1];
			var techErrorsFile = files[2];
			var errorsFile = files[3];
			done = false;
			error = false;
			document.getElementById("tech_error").setAttribute("data-lines", 0);
			document.getElementById("error").setAttribute("data-lines", 0);
			document.getElementById("result").setAttribute("data-lines", 0);
			document.getElementById("detail").setAttribute("data-lines", 0);
			document.getElementById("update").setAttribute("data-lines", 0);
			interval = setInterval(function() {
				getAllOutput(resultsFile, detailsFile, errorsFile, techErrorsFile);
			}, 1000);
		}
	}
	request.open("POST", "hla.php", false);
	request.send(formData);
}

function getOutput(file, divid, divclass, parseLine, headerid) {
	var request = new XMLHttpRequest();
	request.onreadystatechange = function () {
		if (request.readyState == 4 && request.status == 200) {
			var lines = request.responseText.split("\n");
			var	resdiv = document.getElementById(divid);
			var newlines = lines.length-1;
			var oldlines = resdiv.getAttribute("data-lines");
			for (var i=oldlines; i<newlines; ++i) {
				if (divid == "error" || divid == "tech_error") {
					line = lines[i*1];
				} else {
					line = lines[i*1+2];
				}
				if (line.slice(0, 4) == "Done") {
					done = true;
					document.getElementById("status").innerHTML=line;
					return;
				} else if (line) {
					document.getElementById(headerid).className -= "hidden";
					parseLine(line, resdiv);
					if (divid == "tech_error") {
						done = true;
						error = true;
					}
					var lines = resdiv.getAttribute("data-lines")*1;
					resdiv.setAttribute("data-lines", lines+1);
				}
			}
		}
	}
	request.open("GET", file, true);
	request.send(null);
}

function getAllOutput(resultsFile, detailsFile, errorsFile, techErrorsFile) {
	getOutput(techErrorsFile, "tech_error", "tech_error", parseErrorsLine, "tech_error_header");
	getOutput(errorsFile, "error", "error", parseErrorsLine, "error_header");
	if (batchMode) {
		getOutput(resultsFile, "result", "result", parseResultsLine, "result_header");
	} else {
		getOutput(detailsFile, "detail", "detail", parseDetailsLine, "detail_header");
	}
	if (done) {
		if (!error) {
			clearInterval(interval);
			if (batchMode) {
				showDownload(resultsFile);
			}
		} else {
			clearInterval(interval);
			document.getElementById("status").innerHTML="error";	
		}
	}
}

function showDownload(resultsFile) {
	var request = new XMLHttpRequest();
	request.open("POST", "stupid.php", false);
	request.send(null);
	document.getElementById("download").setAttribute("href", resultsFile);
	document.getElementById("download").className -= "hidden";
}

function displayFile (files, textareaid) {
	var reader = new FileReader();
	var textarea = document.getElementById(textareaid);
	textarea.value = "";
	reader.onload = function (f) {
		var text = f.target.result;
		textarea.value += text;
	}
	for (var i=0, f; f=files[i]; i++) {
		reader.readAsText(files[i]);
	}
}

function parseErrorsLine(line, div) {
	var newdiv = document.createElement("div");
	newdiv.className = "error";
	newdiv.appendChild(document.createTextNode(line));
	div.appendChild(newdiv);
}

function parseUpdateLine(line, div) {
	newdiv = document.createElement("div");
	newdiv.innerHTML = line;
	div.appendChild(newdiv);
}

function parseResultsLine (line, resdiv) {
	var div = document.createElement("div");
	div.className = "result";
	var row = line.split(",");

	addField(div, "Sample: ", row[0]);
	addField(div, "Classification: ", row[1]);

	var allAlleles = row[2].split(";");
	var numAlleles = allAlleles.length;
	addField(div, "Possibilities: ", numAlleles);
	addList(div, allAlleles);

	addField(div, "Ambiguous: ", row[3] == 1 ? "yes" : "no");
	addField(div, "Homozygous: ", row[4] == 1 ? "yes" : "no");
	addField(div, "Mismatches: ", row[5]);

	var mismatches = row[6].split(";");
	if (row[6].length > 0) {
		addList(div, mismatches);
	}

	var mislist = new Array();
	for (var i=0; i<mismatches.length; i++) {
		mislist[i] = parseInt(mismatches[i].split(":")[0])-1;
	}
	displaySeq(div, row[7], row[8], row[9], mislist);

	resdiv.appendChild(div);
}

function sumNames(pair) {
	var a0 = pair[0].slice(2).split(":")
	var a1 = pair[1].slice(2).split(":")
	
	while (a0.length < 3) { a0.push("00"); }
	while (a1.length < 3) { a1.push("00"); }

	var sum = a0[0] + a1[0] + a0[1] + a1[1] + a0[2] + a1[2];
	sum = sum.replace(/[NG]/g, "");
	sum = sum.replace(/^0+/g, "");
	return parseInt(sum);
}

function parseDetailsLine(line, resdiv) {
	var div = document.createElement("div");
	div.className = "result";

	var row = line.split(",")

	var mismatches = row[1].length > 0 ? row[1].split(";") : new Array();
	addField(div, "Classification: ", row[0]);
	addField(div, "Mismatches: ", mismatches.length);
	if (row[1].length > 0) {
		addList(div, mismatches);
	}

	var mislist = new Array();
	for (var i=0; i<mismatches.length; i++) {
		mislist[i] = parseInt(mismatches[i].split(":")[0])-1;
	}
	displaySeq(div, row[2], row[3], row[4], mislist);

	var sum = sumNames(row[0].split(" - ")) + mismatches.length*10E15;
	div.setAttribute("data-sum", sum);

	var is_in = false;
	for (var c=resdiv.firstChild; c; c=c.nextSibling) {
		if (c.tagName != "DIV") { continue; }
		if (c.getAttribute("data-sum")*1 > sum) {
			is_in = true;
			resdiv.insertBefore(div, c);
			break;
		}
	}
	if (!is_in) { resdiv.appendChild(div); }
}

function displaySeq(div, exon2, intron, exon3, mismatches) {
	var letter = document.getElementById("letter").value;

	var seqdiv = document.createElement("div");
	seqdiv.className = "sequence";

	var pos = 1;
	exon2_mismatches = mismatches.filter(function (e) {return e < 270;});
	var start = 0;
	for (var i=0; i<exon2_mismatches.length; ++i) {
		end = exon2_mismatches[i];
		pos = addSeq(exon2.slice(start, end), seqdiv, "exon2", pos);
		pos = addSeq(exon2.charAt(end), seqdiv, "exon2 mismatch", pos);
		start = end+1;
	}
	pos = addSeq(exon2.slice(start), seqdiv, "exon2", pos);
	pos = addSeq(intron, seqdiv, "intron", pos);

	exon3_mismatches = mismatches.filter(function(e) {return e >= 270;});
	var start = 0;
	for (i in exon3_mismatches) {
		end = exon3_mismatches[i]-exon2.length;
		if (letter == "A") {
			end -= intron.length;
		}
		pos = addSeq(exon3.slice(start, end), seqdiv, "exon3", pos);
		pos = addSeq(exon3.charAt(end), seqdiv, "exon3 mismatch", pos);
		start = end+1;
	}
	pos = addSeq(exon3.slice(start), seqdiv, "exon3", pos);
	div.appendChild(seqdiv);
}

function addList(div, items) {
	var list = document.createElement("ul");
	for (var i=0; i<items.length; ++i) {
		var item = document.createElement("li");
		item.innerHTML=items[i];
		list.appendChild(item);
	}
	div.appendChild(list);
}

function addField(div, label, value) {
	var newlabel = document.createElement("span");
	newlabel.className = "label";
	newlabel.innerHTML = label;
	div.appendChild(newlabel);
	div.appendChild(document.createTextNode(value));
	div.appendChild(document.createElement("br"));
}

function pad(number, length) {
    var str = '' + number;
    while (str.length < length) {
        str = '0' + str;
    }
    return str;
}

function addSeq(text, div, newclass, pos) {
	var newspan = document.createElement("span");
	newspan.className = newclass;
	for (var i=0; i<text.length; ++i) {
		if (pos % seqWidth == 1) {
			div.appendChild(document.createTextNode(pad(pos, 3) + " "));
		} else if ((pos-1) % seqBlockSize == 0) {
			div.appendChild(newspan);
			div.appendChild(document.createTextNode(" "));
			newspan = document.createElement("span");
			newspan.className = newclass;
		}
		newspan.innerHTML += text.charAt(i);

		if (pos % seqWidth == 0) {
			div.appendChild(newspan);
			div.appendChild(document.createTextNode(" " + pad(pos, 3)));
			div.appendChild(document.createElement("br"));
			newspan = document.createElement("span");
			newspan.className = newclass;
		}

		pos += 1;
	}
	div.appendChild(newspan);
	return pos;
}
