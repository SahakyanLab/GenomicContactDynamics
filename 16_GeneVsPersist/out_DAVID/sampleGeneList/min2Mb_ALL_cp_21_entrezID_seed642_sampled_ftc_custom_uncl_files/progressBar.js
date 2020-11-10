function hideProgress() { 
if (document.getElementById) { // DOM3 = IE5, NS6 
document.getElementById('progress').style.visibility = 'hidden'; 
} 
else { 
if (document.layers) { // Netscape 4 
document.progress.visibility = 'hidden'; 
} 
else { // IE 4 
document.all.progress.style.visibility = 'hidden'; 
} 
} 
} 

function showProgress() { 
if (document.getElementById) { // DOM3 = IE5, NS6 
document.getElementById('progress').style.visibility = 'visible'; 
} 
else { 
if (document.layers) { // Netscape 4 
document.progress.visibility = 'visible'; 
} 
else { // IE 4 
document.all.progress.style.visibility = 'visible'; 
} 
} 
} 
function clearPaste(){}

function setProgress(size) {
 
	if (document.getElementById) { // DOM3 = IE5, NS6 
	document.getElementById('bar').style.width = size*9.8+'px'; 
	document.getElementById('barText').innerHTML = (size)+'%';
	} 
	else { 
		if (document.layers) { // Netscape 4 
		document.bar.width = size*9.8; 
		document.barText.innerHTML = (size)+'%';

		} 
		else { // IE 4 
		document.all.bar.style.width = size*9.8; 
		document.all.barText.innerHTML = (size)+'%';
		} 
	} 
} 

function setProgress2(size) { 
if (document.getElementById) { // DOM3 = IE5, NS6 
document.getElementById('bar').style.width = size*9.8+'px'; 
document.getElementById('barText').innerHTML = size+'%';


} 
else { 
if (document.layers) { // Netscape 4 
document.bar.width = size*2; 
document.barText.innerHTML = size+'%';

} 
else { // IE 4 
document.all.bar.style.width = size*9.8; 
document.all.barText.innerHTML = size+'%';
} 
} 
} 

function setLabel(label) { 
if (document.getElementById) { // DOM3 = IE5, NS6 

document.getElementById('barText').innerHTML = label;


} 
else { 
if (document.layers) { // Netscape 4 
document.barText.innerHTML = label;

} 
else { // IE 4 
document.all.barText.innerHTML = label;
} 
} 
} 


