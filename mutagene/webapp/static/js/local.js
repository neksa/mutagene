async function requestAction(url, options) {
    const response = await fetch(url, options);
    if (!response.ok) {
        const data = await response.json().catch(() => ({}));
        throw new Error(data.error || `Request failed (HTTP ${response.status}).`);
    }
    return response.json();
}

function showRequestError(error) {
    const panel = document.getElementById('requestError');
    panel.textContent = error.message;
    panel.hidden = false;
    panel.scrollIntoView({block: 'nearest'});
}

async function startAnalysis(id) {
    await requestAction(`/api/analyze/${id}`, {method: 'POST'});
    window.location.assign(`/results/${id}`);
}

async function deleteAnalysis(id) {
    if (!window.confirm('Delete this analysis and its files?')) return;
    try {
        await requestAction(`/api/delete/${id}`, {method: 'DELETE'});
        window.location.assign('/history');
    } catch (error) {
        showRequestError(error);
    }
}
