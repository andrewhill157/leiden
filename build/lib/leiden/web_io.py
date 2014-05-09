import urllib


def get_page_html(page_url):
    """
    Returns the html describing the page at the specified URL.

    Args:
        page_url (str): URL to a specified website

    Returns:
        str: HTML describing the specified page

    Raises:
        IOError:
            if an error code is returned
        ValueError:
            if requested URL not found

    """
    try:
        response_code = urllib.urlopen(page_url).code
    except IOError as e:
        raise ValueError('Requested URL not found. Please ensure http:// is included.')

    if response_code >= 400:
        raise IOError('Requested URL could not be reached. Response code: ' + response_code)

    return urllib.urlopen(page_url).read()