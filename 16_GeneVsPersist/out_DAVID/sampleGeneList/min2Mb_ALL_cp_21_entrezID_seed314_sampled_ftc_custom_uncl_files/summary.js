// This is variable for storing callback function
var ae_cb = null;
 

// this is a simple function-shortcut
// to avoid using lengthy document.getElementById
function ae$(a) { return document.getElementById(a); }
 
// This is a main ae_prompt function
// it saves function callback 
// and sets up dialog
function ae_prompt(cb, q, a) {

	ae_cb = cb;
	//ae$('aep_t').innerHTML = document.domain + ' prompt:';
	ae$('aep_prompt').innerHTML = q;
	ae$('aep_text').value = a;
	ae$('aep_ovrl').style.display = ae$('aep_ww').style.display = '';
	ae$('aep_text').focus();
	ae$('aep_text').select();

}
 
// This function is called when user presses OK(m=0) or Cancel(m=1) button
// in the dialog. You should not call this function directly.
function ae_clk(m) {
	// hide dialog layers 
	ae$('aep_ovrl').style.display = ae$('aep_ww').style.display = 'none';
	if (!m)  
		ae_cb(null);  // user pressed cancel, call callback with null
	else
		ae_cb(ae$('aep_text').value); // user pressed OK 
}



function checkDownloadForm()
{
	var FM=document.getElementById("frmDBDownload");
	
	for (var i=0;i<FM.termFileList.length;i++)
	{
		if(FM.termFileList[i].checked)
		{
			FM.submit();
			return;
		}

	}

	alert("You must select at least one file to download.");
	return;

}
function checkEmail()
{
	var FM=document.getElementById("frmDownloadLogin");
	if(FM.email==null||FM.email.value=="")
	{
		alert("You must enter a valid email address.");
		return;
	}
	FM.submit();
}

function register()
{
	var FM=document.getElementById("frmDBDownload");

	if(FM.first==null||FM.first.value=="")
	{
		alert("Please enter your first name.");
		return;
	}
	if(FM.last==null||FM.last.value=="")
	{
		alert("Please enter your last name.");
		return;
	}
	if(FM.org==null||FM.org.value=="")
	{
		alert("Please enter the name of your organization.");
		return;
	}
	if(FM.email==null||FM.email.value=="")
	{
		alert("You must enter a valid email address.");
		return;
	}
var emailStr=FM.email.value.toUpperCase();
if(emailStr.indexOf("HOTMAIL")>0||emailStr.indexOf("GMAIL.COM")>0||emailStr.indexOf("YAHOO")>0)
	{
		alert("You must enter an email address from your academic, government, or commercial organization.\nEmail addresses from public providers will not be accepted.");
		return;
	}


var checked=false;
for(var i = 0; i < FM.entity.length; i++) 
{
		if(FM.entity[i].checked) 
		{
			checked=true;
		}
}

if(!checked)
	{
		alert("Please choose the type of organization to which you belong.");
		return;
	}



FM.action='downloadRegistration.jsp';
FM.submit();
}

function submitConvertedList()
{
	//alert(document.getElementById("VTO").value);
	//alert(document.getElementById("convertedIds").value);
	
var win=false;
try
   {
   	win=window.opener==null;
	if(!win)
	{
		win=window.opener.closed;
	}

   }
catch(err)
   {
   	win=true;
   }


//alert(win);


	if(win)
	{
		alert("You have closed the main DAVID window.  You will now be redirected to start the conversion process again.");
		var FM=document.getElementById("FRM");
		FM.action='conversion.jsp?VFROM='+document.getElementById("VFROM").value;FM.target='';FM.submit();
	}
	else
	{
	if(document.getElementById("convertedIds").value==null||document.getElementById("convertedIds").value=="")
	{
		alert("You must first convert your ids.");
		return;
	}
	var FM = window.opener.document.forms[0];
	FM.ids.value=document.getElementById("convertedIds").value;
	FM.idType.value=document.getElementById("VTO").value;
       //alert(FM.idType.value);
	FM.uploadType.value="list";

	FM.Mode.value="paste";
	
	FM.action="tools.jsp";
	FM.target="";
	FM.submit();
	window.opener.focus();
	}
	//alert(FM.idType.value);



}
function submitConvertedPopulation()
{
	//alert(document.getElementById("VTO").value);
	//alert(document.getElementById("convertedIds").value);
	
var win=false;
try
   {
   	win=window.opener==null;
	if(!win)
	{
		win=window.opener.closed;
	}

   }
catch(err)
   {
   	win=true;
   }


//alert(win);


	if(win)
	{
		alert("You have closed the main DAVID window.  You will now be redirected to start the conversion process again.");
		var FM=document.getElementById("FRM");
		FM.action='conversion.jsp?VFROM='+document.getElementById("VFROM").value;FM.target='';FM.submit();
	}
	else
	{
	if(document.getElementById("convertedIds").value==null||document.getElementById("convertedIds").value=="")
	{
		alert("You must first convert your ids.");
		return;
	}
	var FM = window.opener.document.forms[0];
	FM.ids.value=document.getElementById("convertedIds").value;
	FM.idType.value=document.getElementById("VTO").value;
       //alert(FM.idType.value);
	FM.uploadType.value="population";

	FM.Mode.value="paste";
	
	FM.action="tools.jsp";
	FM.target="";
	FM.submit();
	window.opener.focus();
	}
	//alert(FM.idType.value);



}

function checkDefaults()
{
	var state=document.forms['summaryForm'].elements['checkDefault'].checked;
	
	var defaultsArray=defaults.split(",");
	if(!document.forms['summaryForm'])
	  {
		return;
	  }	
	  //alert(defaults);
	  for (i=0; i<defaultsArray.length; i++) 
	  {
		  if(defaultsArray[i]!='')
		  {
			document.forms['summaryForm'].summaryList[defaultsArray[i]].checked=state;
		  }
		
	  }

	countAll();
}
function uncheckAll()
{

	  if(!document.forms['summaryForm'])
	  {
		return;
	  }	

	  for (i=0; i<document.forms['summaryForm'].summaryList.length; i++) 
	  {
		document.forms['summaryForm'].summaryList[i].checked=false;
		
	  }


	
	document.forms['summaryForm'].elements['checkDefault'].checked=false;
	countAll();
}

function submitConvertedGene(){
				ae_prompt(submitConvertedGeneAction,'Enter a new name for this list','new_converted_list');

				
			}
function submitConvertedGeneAction(n)
{
	if(n==null)return;
var listname =n;
				
				var win=false;
				try{
				   	win=window.opener==null;
					if(!win){
						win=window.opener.closed;
					}
				
				}catch(err){
				   	win=true;
				}
				if(win){
						alert("You have closed the main DAVID window.\nThis will redirect you to DAVID window and start the conversion process again.");
						window.location = 'convesion.jsp';
				}else{
					var FM = window.opener.document.forms[0];
					FM.convertedListName.value = listname;	
					FM.action="tools.jsp";
					FM.target="";
					FM.submit();
					window.opener.focus();
				}

	
}

			
			function submitConvertedBg(){
				ae_prompt(submitConvertedBgAction,'Enter a new name for this background list','new_converted_bg');
				
			}
		
function submitConvertedBgAction(n)
{
	if(n==null)return;
var listname =n;
				
				var win=false;
				try{
				   	win=window.opener==null;
					if(!win){
						win=window.opener.closed;
					}
				
				}catch(err){
				   	win=true;
				}
				if(win){
						alert("You have closed the main DAVID window.\nThis will redirect you to DAVID window and start the conversion process again.");
						window.location = 'conversion.jsp';
				}else{
					var FM = window.opener.document.forms[0];
					FM.convertedPopName.value = listname;	
					FM.action="tools.jsp";
					FM.target="";
					FM.submit();
					window.opener.focus();
				}	

	
}
function createSublist()
{
	
ae_prompt(createSublistAction,'Enter a name for this sublist:','sublist'); 

 
}

function createSublistAction(n)
{
	if(n==null)return;
	var  sublistStr='';
	var form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);
	}
	var x=0; 
	var boxes=document.forms[form].elements['CKDL'];

	if(!boxes)
	{
		alert('no checkboxes');
		return;
	}
	var count=0;
	for (var i=0; i<boxes.length; i++) 
	{ 
		if(boxes[i].checked)
		{
			sublistStr=sublistStr+','+boxes[i].value;
			count++;
		}
			
	}
	if(count==0)
	{
		alert('no boxes checked');
		return;
	}
	//alert(sublistStr);
	//alert(form);


if(form==0)
{


	window.opener.document.forms[0].rowids.value=sublistStr;
	window.opener.document.forms[0].sublist.value=n;
	window.opener.document.forms[0].submit();
	window.opener.focus();

}
else
{
	document.forms[0].rowids.value=sublistStr;
	document.forms[0].sublist.value=n;
	document.forms[0].submit();
}

	
}
function geneReport(dids)
{
	var form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);
	}
	//alert(document.forms[form].id);
	setRowids(dids);
	setAction('geneReport.jsp');

	//alert(document.getElementById("rowids").value);

	document.forms[form].submit();
}
function geneReportFull(dids)
{
	var form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);
	}
	//alert(document.forms[form].id);
	setRowids(dids);
	setAction('geneReportFull.jsp');

	//alert(document.getElementById("rowids").value);

	document.forms[form].submit();



}
function termReport(ids,didCount)
{
	
	var form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);
	}

	setRowids(ids);setAction('termReport.jsp?didCount='+didCount);document.forms[form].submit();

}
function termReport2()
{
	
	var form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);
	}

	setAction('termReport2.jsp');document.forms[form].submit();

}


function submitForm(action,target,form)
{
	form.action=action;
	form.target=target;
	form.submit();

}
function submitRelatedGenes(action,target,form)
{
	//if(form.list.value==-1)
	//{
	//	alert("This is an intensive operation which may take several minutes.");
	//}
	form.action=action;
	form.target=target;
	form.submit();

}
function gene2gene()
{
	var wnd = window.open("gene2gene.jsp?stringency="+document.forms['gene2gene'].stringency.value+"&kappa="+document.forms['gene2gene'].kappa.value+"&overlap="+document.forms['gene2gene'].overlap.value+"&initialSeed="+document.forms['gene2gene'].initialSeed.value+"&frequency="+document.forms['gene2gene'].frequency.value+"&finalSeed="+document.forms['gene2gene'].finalSeed.value+"&linkage="+document.forms['gene2gene'].linkage.value, "_self", "width=700, height=500, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();

}
function setCheckCount(grouptd,count){
					
			
			var obj=document.getElementById(grouptd);
			//alert(obj.innerHTML);
			obj.innerHTML="&nbsp;("+count+" selected)";
		}
function countGroup(GroupNode) {
		//alert(GroupNode.parentNode.parentNode.parentNode.parentNode.id+"->"+GroupNode.parentNode.parentNode.parentNode.id+"->"+GroupNode.parentNode.parentNode.id+"->"+GroupNode.parentNode.id+"->"+GroupNode.id);
		if(!document.forms['summaryForm'])
		{
		return;
		}
	var CB = document.forms['summaryForm'].elements[GroupNode];
	if(CB)
	{
					    var count=0;
				
				for(var i = 0; i < CB.length; i++)
			{
				if(CB[i].checked)
				{
					count++;
				}
			}
			 setCheckCount(GroupNode+"td",count);
	}
}
function countAll()
{
	
          var group='';
	  if(!document.forms['summaryForm'])
	  {
		return;
	  }	

	  for (i=0; i<document.forms['summaryForm'].summaryList.length; i++) 
	  {
		if(group!=document.forms['summaryForm'].summaryList[i].id)
		{
			group=document.forms['summaryForm'].summaryList[i].id;
			//alert(group);
			countGroup(group);
		}
		
	  }

		
		
}
function getCurrentList()
{

		
		if(!document.forms['summaryForm'])
		{
			return;
		}
		return document.getElementById("currentList").value;

}
function getTerms()
{
	var  termStr='';
	var selectCount=0;

	  if(!document.forms['summaryForm'])
	  {
		return;
	  }	
	  
	  for (i=0; i<document.forms['summaryForm'].summaryList.length; i++) 
	  {
	    	if(document.forms['summaryForm'].summaryList[i].checked)
		{
			if(termStr!='')
			{
				termStr=termStr+','+document.forms['summaryForm'].summaryList[i].value;
			}
			else
			{
				termStr=document.forms['summaryForm'].summaryList[i].value;
			
			}
			selectCount++;
		}
	  }

	
	if(selectCount>0)
	{		 
		//alert(termStr);
		return termStr;
	}
	else
	{
		alert('You must select an annotation category.');
		return null;
	}
			 
}


function Details(group, term, idt){

	var wnd = window.open("geneReport.jsp?TX='9606:Homo sapiens'&RF=Y&LS=demolist1&GR=" + group + "&TM=" + term + "&ID=" + idt, "", "width=700, height=500, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}
function twoDView(index,type){

	var wnd = window.open("2DView.jsp?majority=true&index=" + index+"&clusterType=" + type, "", "");
	wnd.focus();
}
function setAction(action)
{
	
	form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);

	}

	document.forms[form].action=action;
	//alert(action);
	//alert(document.forms[form].id);

	if(action!="term2term.jsp"&&action!="gene2gene.jsp"&&action!="relatedTerms.jsp"&&action!="relatedGenes.jsp")
	{
		document.forms[form].target="_blank";
	}
	else
	{
		document.forms[form].target="_self";
	}

}
function setRowids(rowids)
{
	var form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);
	}
	
	//document.getElementById("rowids").value=rowids;
document.forms[form].rowids.value=rowids;

	//alert(document.getElementById("rowids").value);
	//alert(document.forms[form].id);

}
function setHeading(headingStr)
{
	var form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);
	}

	document.forms[form].heading.value=headingStr;
}
function setClusterSort(sort)
{
	var form=0;
	if(document.forms.length>1)
	{
		form=(document.forms.length-1);
	}

	document.forms[0].clusterSort.value=sort;
}

function RT(cat, trm, UT){
if (cat== "GOTERM_BP_ALL") cat="1";
else if (cat== "GOTERM_MF_ALL") cat="2";
else if (cat== "GOTERM_CC_ALL") cat="3";
else if (cat== "KEGG_PATHWAY") cat="4";
else if (cat== "BIOCARTA") cat="5";
else if (cat== "SP_PIR_KEYWORDS") cat="6";
else if (cat== "BBID") cat="7";
else if (cat== "SMART_NAME") cat="8";
else if (cat== "GENETIC_ASSOCIATION_DB") cat="9";
else if (cat== "UP_SEQ_FEATURE") cat="10";
else if (cat== "COG_KOG_ONTOLOGY") cat="11";
else if (cat== "OMIM_PHENOTYPE") cat="12";
else if (cat== "INTERPRO_NAME") cat="13";
else if (cat== "PIR_SUPERFAMILY_NAME") cat="14";
else cat="0";
	if(cat == "0"){
		alert("The category is not supported for related Search tool.");
		return;
	}
	var wnd = window.open("relatedTerms.jsp?LT="+UT+"&CAT=" + cat + "&TRM=" + trm, "RT", "width=880, height=550, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}
function relatedTerms(termId){
	var wnd = window.open("relatedTerms.jsp?id=" + termId, "Related Terms", "width=880, height=550, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();

}
function relatedGenes(rowId){
	var wnd = window.open("relatedGenes.jsp?id=" + rowId, "Related Genes", "width=880, height=550, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();

}


function LinkOut(cat, trm){		
	var url = "";
	
	if(cat.toUpperCase() == "GOTERM_BP_ALL"){
		url = "http://www.ebi.ac.uk/ego/GSearch?query=" + trm+ "&mode=name";
	}
	if(cat.toUpperCase() == "UP_SEQ_FEATURE"){
		url = "http://www.pir.uniprot.org/";
	}
	if(cat.toUpperCase() == "GOTERM_MF_ALL"){
		url = "http://www.ebi.ac.uk/ego/GSearch?query=" + trm + "&mode=name";
	}
	
	if(cat.toUpperCase() == "SP_PIR_KEYWORDS"){
		url = "http://us.expasy.org/cgi-bin/get-entries?KW=" + trm;
	}	 
	
	if(cat.toUpperCase() == "BIOCARTA"){
			url = "biocarta.jsp?PATH=" + trm.replace(":", "$") + "&TB=demolist1";
	}	
	
	if(cat.toUpperCase() == "KEGG_PATHWAY"){
			url = "kegg.jsp?PATH=" + trm.replace(":", "$") + "&TB=demolist1";
	}	
	
	
	if(cat.toUpperCase() == "GOTERM_CC_ALL"){
				url = "http://www.ebi.ac.uk/ego/GSearch?query=" + trm + "&mode=name";
	}
	
	if(cat.toUpperCase() == "BBID"){
					url = "http://bbid.grc.nia.nih.gov/geneimages/" + trm + ".jpeg";
	}
	
	if(cat.toUpperCase() == "PIR_SUPERFAMILY_NAME"){
					url = "http://pir.georgetown.edu/cgi-bin/ipcSF?id=" + trm.substring(0, trm.indexOf(":"));
	}
	
	if(cat.toUpperCase() == "SMART_NAME"){
					url = "http://smart.embl-heidelberg.de/smart/do_annotation.pl?BLAST=DUMMY&ACC=" + trm.substring(0, trm.indexOf(":"));
	}
	
	if(cat.toUpperCase() == "INTERPRO_NAME"){
					url = "http://www.ebi.ac.uk/interpro/DisplayIproEntry?ac=" + trm.substring(0, trm.indexOf(":")) ;
	}
	
	if(cat.toUpperCase() == "SP_PIR_KEYWORDS"){
					url = "http://us.expasy.org/cgi-bin/get-entries?KW=" + trm;
	}
	
	if(cat.toUpperCase() == "OMIM_PHENOTYPE"){
					url = "http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=" + trm.substring(0, trm.indexOf(":"));
	}
	
	if(cat.toUpperCase() == "COG_KOG_ONTOLOGY"){
					url = "http://www.ncbi.nlm.nih.gov/COG/new/ ";
	}
	if(cat.toUpperCase() == "GENETIC_ASSOCIATION_DB"){
					url = "http://geneticassociationdb.nih.gov/";
	}	
	
	
	if(url != ""){
		var wnd = window.open(url, "LKOT", "width=780, height=600, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
		wnd.focus();
	}
		
}
function CreateNRChart(){
	
		if(getTerms()!=null)
		{
			var wnd = window.open("term2term.jsp?annot=" + getTerms()+"&currentList=" + getCurrentList(), "", "width=980, height=620, toolbar=yes, menubar=yes, resizable=yes, scrollbars=yes");
			wnd.focus();
		}

		

		return;

}
function CreateCluster(){
	
		if(getTerms()!=null)
		{
			var wnd = window.open("cluster.jsp?annot=" + getTerms()+"&currentList=" + getCurrentList(), "", "width=980, height=620, toolbar=yes, menubar=yes, resizable=yes, scrollbars=yes");
			wnd.focus();
		}
		return;

}
function CreateTable()
{

	
	if(getTerms()!=null)

	{

		var wnd = window.open("annotationReport.jsp?annot=" + getTerms()+"&currentList=" + getCurrentList(), "", "width=800, height=620, toolbar=yes, menubar=yes, resizable=yes, scrollbars=yes");
		wnd.focus();
	}
		

		return;
}
function knowledgebaseReport()
{

	
	if(getTerms()!=null)

	{

		var wnd = window.open("knowledgebaseReport.jsp?annot=" + getTerms()+"&currentList=" + getCurrentList(), "", "width=800, height=620, toolbar=yes, menubar=yes, resizable=yes, scrollbars=yes");
		wnd.focus();
	}
		

		return;
}
function exportTable()
{
	
	if(getTerms()!=null)
	{

		var wnd = window.open("annotationExport.jsp?annot=" + getTerms(), "", "width=800, height=620, toolbar=yes, menubar=yes, resizable=yes, scrollbars=yes");
		wnd.focus();
	}
		

		return;
}

function CreateChart()
{
	
	if(getTerms()!=null)
	{

		var wnd = window.open("chartReport.jsp?annot=" + getTerms()+"&currentList=" + getCurrentList(), "", "width=800, height=620, toolbar=yes, menubar=yes, resizable=yes, scrollbars=yes");
		wnd.focus();
	}
		

		return;


}
function ShowGeneList()
{
	
	window.location = "list.jsp";
		


}
function exportChart()
{
	
	if(getTerms()!=null)
	{

		var wnd = window.open("chartExport.jsp?annot=" + getTerms(), "", "width=800, height=620, toolbar=yes, menubar=yes, resizable=yes, scrollbars=yes");
		wnd.focus();
	}

		return;


}


var plus = new Image();
  plus.src = "images/plus.gif";
var minus = new Image();
  minus.src = "images/minus.gif";
var item = new Image();
  item.src = "images/item.gif";

function door(obj){

  obj = obj.parentNode;
  obj.className = (obj.className == 'open') ? 'closed' : 'open';
  return false;
}

function getAnnotationReport(annot)
{
	var wnd = window.open("annotationReport.jsp?annot=" + annot, "", "width=720, height=350, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}
function getKnowledgebaseReport(annot)
{
	var wnd = window.open("knowledgebaseReport.jsp?annot=" + annot, "", "width=720, height=350, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}
function getChartReport(annot)
{
	var wnd = window.open("chartReport.jsp?annot=" + annot, "", "width=720, height=350, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}

function OpenHelp(part){
	var wnd = window.open("helps/linear_search.html#" + part, "Details", "width=720, height=350, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}

function OpenHelp2(part){
	var wnd = window.open("helps/functional_annotation.html#" + part, "Details", "width=720, height=350, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}
function CombinedReport(usertable,did){
	var wnd = window.open("geneReportFull.jsp?LT=" + usertable + "&ID=" + did, "Report", "width=640, height=500, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}


function RG(val){	
	var wnd = window.open("relatedGenes.jsp?LT=demolist1&ID=" + val, "", "width=880, height=360, toolbar=no, menubar=no, resizable=yes, scrollbars=yes");
	wnd.focus();
}
