;(function($) {
	$.loadBinary = function(url, callback) {
		var data = null;
		var xhr = new XMLHttpRequest();
		var async = callback != null;
		xhr.open('GET', url, async);
		xhr.responseType = 'arraybuffer';
		if (xhr.overrideMimeType) {
			xhr.overrideMimeType('text/plain; charset=x-user-defined');
		}
		xhr.onload = function(e) {
			var buff = xhr.response;
			data = new Uint8Array(buff);
			if (callback != null) {
				callback(data);
			}
		};
		xhr.send();
		return data;
	};
})($);
