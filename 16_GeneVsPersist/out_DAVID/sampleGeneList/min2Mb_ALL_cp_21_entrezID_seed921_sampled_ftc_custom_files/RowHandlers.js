// Adds onmouseover, onmouseout, and onclick handlers to each table row.  The onmouseover handler changes the row's class attribute to
// rowMouseOver.  The onmouseout handler changes it back.  The onclick function makes a request for the specified url, including the
// innerHTML of the specified column as a request parameter.
function addRowHandlers(tableId, rowClassName, url, paramName, columnIndex) {
    var previousClass = null;
    var table = document.getElementById(tableId);
    var rows = table.getElementsByTagName("tr");
    for (i = 1; i < rows.length; i++) {
        rows[i].onmouseover = function () {
            previousClass = this.className;
            this.className = this.className + " " + rowClassName ;
        };
        
        rows[i].onmouseout = function () {
            this.className = previousClass;
        };
        
       // rows[i].onclick = function () {
       //     var cell = this.getElementsByTagName("td")[columnIndex];
       //     var paramValue = cell.innerHTML;
       //     location.href = url + "?" + paramName + "=" + paramValue;
       // };
    }
}

